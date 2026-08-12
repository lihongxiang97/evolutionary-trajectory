#!/usr/bin/env python3
"""Quantify concordance and age-effect replication between pear haplotypes."""
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
from statsmodels.genmod.families import Binomial


TYPES = ["TD", "PD", "TRD"]


def as_bool(series):
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def bh(values):
    values = np.asarray(values, dtype=float)
    result = np.full_like(values, np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    if not len(valid):
        return result
    order = valid[np.argsort(values[valid])]
    ranked = values[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(order) / np.arange(1, len(order) + 1))[::-1]
    )[::-1]
    result[order] = np.minimum(adjusted, 1)
    return result


def cohen_kappa(left, right):
    left = np.asarray(left, dtype=bool)
    right = np.asarray(right, dtype=bool)
    observed = np.mean(left == right)
    p_left = np.mean(left)
    p_right = np.mean(right)
    expected = p_left * p_right + (1 - p_left) * (1 - p_right)
    return (observed - expected) / (1 - expected) if expected < 1 else np.nan


def odds_ratio_ci(table):
    array = np.asarray(table, dtype=float)
    corrected = array + 0.5 if np.any(array == 0) else array
    log_or = math.log(
        corrected[0, 0]
        * corrected[1, 1]
        / (corrected[0, 1] * corrected[1, 0])
    )
    standard_error = math.sqrt(np.sum(1 / corrected))
    return (
        math.exp(log_or),
        math.exp(log_or - 1.96 * standard_error),
        math.exp(log_or + 1.96 * standard_error),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-features", required=True, type=Path)
    parser.add_argument("--age-features", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    all_features = pd.read_csv(args.all_features, sep="\t")
    age = pd.read_csv(args.age_features, sep="\t")

    correlations = []
    for metric, left, right in [
        ("fruit_expression", "hapA_TPM", "hapB_TPM"),
        ("compartment_PC1", "hapA_PC1", "hapB_PC1"),
    ]:
        block = all_features[[left, right]].replace([np.inf, -np.inf], np.nan).dropna()
        if metric == "fruit_expression":
            x = np.log2(block[left] + 1)
            y = np.log2(block[right] + 1)
        else:
            x, y = block[left], block[right]
        result = spearmanr(x, y)
        correlations.append(
            {
                "metric": metric,
                "n": len(block),
                "spearman_rho": result.statistic,
                "pvalue": result.pvalue,
            }
        )
    for metric, left, right in [
        ("TSS_ATAC", "hapA_mean_ATAC_TSS_2kb", "hapB_mean_ATAC_TSS_2kb"),
        (
            "marginal_contact",
            "hapA_log2_contact_CPM_plus_1",
            "hapB_log2_contact_CPM_plus_1",
        ),
    ]:
        block = age[[left, right]].replace([np.inf, -np.inf], np.nan).dropna()
        result = spearmanr(block[left], block[right])
        correlations.append(
            {
                "metric": metric,
                "n": len(block),
                "spearman_rho": result.statistic,
                "pvalue": result.pvalue,
            }
        )
    pd.DataFrame(correlations).to_csv(
        args.outdir / "continuous_feature_correlations.tsv", sep="\t", index=False
    )

    concordance = []
    for metric, left, right, transform in [
        (
            "A_compartment",
            "hapA_compartment",
            "hapB_compartment",
            lambda frame, column: frame[column].eq("A"),
        ),
        (
            "any_loop",
            "hapA_any_loop",
            "hapB_any_loop",
            lambda frame, column: as_bool(frame[column]),
        ),
        (
            "GGL",
            "hapA_GGL",
            "hapB_GGL",
            lambda frame, column: as_bool(frame[column]),
        ),
        (
            "DGL",
            "hapA_DGL",
            "hapB_DGL",
            lambda frame, column: as_bool(frame[column]),
        ),
    ]:
        block = all_features[[left, right]].dropna()
        x, y = transform(block, left), transform(block, right)
        concordance.append(
            {
                "metric": metric,
                "n": len(block),
                "agreement": np.mean(x.to_numpy() == y.to_numpy()),
                "cohen_kappa": cohen_kappa(x, y),
                "hapA_positive_fraction": np.mean(x),
                "hapB_positive_fraction": np.mean(y),
            }
        )
    block = age[
        age["hapA_TAD_location"].isin(["Boundary", "Interior"])
        & age["hapB_TAD_location"].isin(["Boundary", "Interior"])
    ]
    x = block["hapA_TAD_location"].eq("Boundary")
    y = block["hapB_TAD_location"].eq("Boundary")
    concordance.append(
        {
            "metric": "TAD_boundary",
            "n": len(block),
            "agreement": np.mean(x.to_numpy() == y.to_numpy()),
            "cohen_kappa": cohen_kappa(x, y),
            "hapA_positive_fraction": np.mean(x),
            "hapB_positive_fraction": np.mean(y),
        }
    )
    pd.DataFrame(concordance).to_csv(
        args.outdir / "binary_feature_concordance.tsv", sep="\t", index=False
    )

    effects = []
    continuous = {
        "expression": ("log2_TPM_plus_1", False),
        "ATAC": ("mean_ATAC_TSS_2kb", False),
        "PC1": ("PC1", False),
        "contact": ("log2_contact_CPM_plus_1", False),
    }
    for haplotype in ("hapA", "hapB"):
        for duplication_type in TYPES:
            subset = age[age["duplication_type"].eq(duplication_type)]
            for metric, (suffix, _) in continuous.items():
                column = f"{haplotype}_{suffix}"
                young = subset.loc[subset["age"].eq("young"), column].dropna()
                old = subset.loc[subset["age"].eq("old"), column].dropna()
                result = mannwhitneyu(young, old, alternative="two-sided")
                rank_biserial = 2 * result.statistic / (len(young) * len(old)) - 1
                effects.append(
                    {
                        "haplotype": haplotype,
                        "duplication_type": duplication_type,
                        "metric": metric,
                        "effect_type": "rank_biserial",
                        "effect": rank_biserial,
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "n_young": len(young),
                        "n_old": len(old),
                        "median_young": young.median(),
                        "median_old": old.median(),
                        "pvalue": result.pvalue,
                    }
                )

            binary_specs = [
                (
                    "A_compartment",
                    subset[f"{haplotype}_compartment"].eq("A"),
                    subset[f"{haplotype}_compartment"].isin(["A", "B"]),
                ),
                (
                    "TAD_boundary",
                    subset[f"{haplotype}_TAD_location"].eq("Boundary"),
                    subset[f"{haplotype}_TAD_location"].isin(["Boundary", "Interior"]),
                ),
                (
                    "any_loop",
                    as_bool(subset[f"{haplotype}_any_loop"]),
                    subset[f"{haplotype}_any_loop"].notna(),
                ),
            ]
            for metric, positive, eligible in binary_specs:
                labels = subset.loc[eligible, "age"]
                outcomes = positive[eligible]
                young_positive = int(outcomes[labels.eq("young")].sum())
                young_negative = int(labels.eq("young").sum() - young_positive)
                old_positive = int(outcomes[labels.eq("old")].sum())
                old_negative = int(labels.eq("old").sum() - old_positive)
                table = [
                    [young_positive, young_negative],
                    [old_positive, old_negative],
                ]
                _, pvalue = fisher_exact(table)
                odds_ratio, ci_low, ci_high = odds_ratio_ci(table)
                effects.append(
                    {
                        "haplotype": haplotype,
                        "duplication_type": duplication_type,
                        "metric": metric,
                        "effect_type": "odds_ratio",
                        "effect": odds_ratio,
                        "ci_low": ci_low,
                        "ci_high": ci_high,
                        "n_young": young_positive + young_negative,
                        "n_old": old_positive + old_negative,
                        "median_young": np.nan,
                        "median_old": np.nan,
                        "pvalue": pvalue,
                    }
                )
    effects = pd.DataFrame(effects)
    effects["p_BH"] = bh(effects["pvalue"])
    effects.to_csv(
        args.outdir / "haplotype_age_effects.tsv", sep="\t", index=False
    )

    pivot = effects.pivot_table(
        index=["duplication_type", "metric", "effect_type"],
        columns="haplotype",
        values="effect",
    ).reset_index()
    pivot["hapA_direction"] = np.where(
        pivot["effect_type"].eq("odds_ratio"),
        np.sign(np.log(pivot["hapA"])),
        np.sign(pivot["hapA"]),
    )
    pivot["hapB_direction"] = np.where(
        pivot["effect_type"].eq("odds_ratio"),
        np.sign(np.log(pivot["hapB"])),
        np.sign(pivot["hapB"]),
    )
    pivot["direction_concordant"] = pivot["hapA_direction"].eq(
        pivot["hapB_direction"]
    )
    pivot.to_csv(
        args.outdir / "age_effect_direction_concordance.tsv",
        sep="\t",
        index=False,
    )

    adjusted = []
    frame = age.copy()
    frame["age_young"] = frame["age"].eq("young").astype(int)
    frame["log_gene_length"] = np.log1p(frame["hapB_gene_length"])
    frame["log_local_density"] = np.log1p(frame["hapB_local_gene_density_200kb"])
    frame["A_compartment"] = frame["hapB_compartment"].eq("A").astype(int)
    frame["TAD_boundary"] = frame["hapB_TAD_location"].eq("Boundary").astype(int)
    frame["any_loop_binary"] = as_bool(frame["hapB_any_loop"]).astype(int)
    base = (
        "age_young + C(hapB_chrom) + log_gene_length + hapB_gene_body_GC + "
        "hapB_exon_count + log_local_density"
    )
    for duplication_type in TYPES:
        subset = frame[frame["duplication_type"].eq(duplication_type)].copy()
        for metric, outcome, include_expression, family in [
            ("expression", "hapB_log2_TPM_plus_1", False, "ols"),
            ("ATAC", "hapB_mean_ATAC_TSS_2kb", False, "ols"),
            ("PC1", "hapB_PC1", True, "ols"),
            ("contact", "hapB_log2_contact_CPM_plus_1", True, "ols"),
            ("A_compartment", "A_compartment", True, "binomial"),
            ("TAD_boundary", "TAD_boundary", True, "binomial"),
            ("any_loop", "any_loop_binary", True, "binomial"),
        ]:
            columns = [
                outcome,
                "age_young",
                "hapB_chrom",
                "log_gene_length",
                "hapB_gene_body_GC",
                "hapB_exon_count",
                "log_local_density",
                "hapB_log2_TPM_plus_1",
            ]
            # The expression outcome is also the optional expression covariate.
            # De-duplicate the column list so patsy receives one numeric response
            # rather than two identically named response columns.
            columns = list(dict.fromkeys(columns))
            model_data = subset[columns].replace([np.inf, -np.inf], np.nan).dropna()
            formula = f"{outcome} ~ {base}"
            if include_expression and outcome != "hapB_log2_TPM_plus_1":
                formula += " + hapB_log2_TPM_plus_1"
            try:
                if family == "ols":
                    model = smf.ols(formula, data=model_data).fit(cov_type="HC3")
                    effect = model.params["age_young"]
                    low, high = model.conf_int().loc["age_young"]
                else:
                    model = smf.glm(
                        formula, data=model_data, family=Binomial()
                    ).fit(cov_type="HC3")
                    effect = math.exp(model.params["age_young"])
                    interval = model.conf_int().loc["age_young"]
                    low, high = math.exp(interval.iloc[0]), math.exp(interval.iloc[1])
                adjusted.append(
                    {
                        "duplication_type": duplication_type,
                        "metric": metric,
                        "model": family,
                        "effect": effect,
                        "ci_low": low,
                        "ci_high": high,
                        "pvalue": model.pvalues["age_young"],
                        "n": int(model.nobs),
                    }
                )
            except Exception as error:
                adjusted.append(
                    {
                        "duplication_type": duplication_type,
                        "metric": metric,
                        "model": family,
                        "effect": np.nan,
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "pvalue": np.nan,
                        "n": len(model_data),
                        "error": repr(error),
                    }
                )
    adjusted = pd.DataFrame(adjusted)
    adjusted["p_BH"] = bh(adjusted["pvalue"])
    adjusted.to_csv(
        args.outdir / "hapB_adjusted_age_models.tsv", sep="\t", index=False
    )

    with (args.outdir / "cross_haplotype_robustness_summary.tsv").open("w") as output:
        output.write("metric\tvalue\n")
        output.write(f"core_RBH_genes\t{len(all_features)}\n")
        output.write(f"age_classified_RBH_genes\t{len(age)}\n")
        output.write(
            f"age_effect_direction_concordance\t"
            f"{pivot['direction_concordant'].mean():.6f}\n"
        )
        output.write(
            f"age_effect_direction_concordant_tests\t"
            f"{int(pivot['direction_concordant'].sum())}\n"
        )
        output.write(f"age_effect_tests\t{len(pivot)}\n")


if __name__ == "__main__":
    main()
