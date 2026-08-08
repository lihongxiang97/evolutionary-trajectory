#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TRAJ = ROOT / "08.trajectories/results/trajectory_classifications.tsv"
KAKS_DIR = ROOT / "06.kaks/SSD/output"
OUT = Path(__file__).resolve().parent / "source"
OUT.mkdir(parents=True, exist_ok=True)

traj = pd.read_csv(TRAJ, sep="\t")
traj = traj[
    (traj["analysis_status"] == "classified")
    & traj["two_method_agreement"].isin(
        [
            "Conservation",
            "Neofunctionalization(Child)",
            "Neofunctionalization(Parent)",
            "Specialization",
        ]
    )
].copy()

rows = []
for path in KAKS_DIR.glob("*.kaks"):
    try:
        table = pd.read_csv(path, sep="\t")
    except Exception:
        continue
    if table.empty:
        continue
    row = table.iloc[0]
    rows.append(
        {
            "pair": str(row.get("Sequence", path.name.split(".cds_aln")[0])),
            "method": row.get("Method"),
            "Ka": pd.to_numeric(row.get("Ka"), errors="coerce"),
            "Ks": pd.to_numeric(row.get("Ks"), errors="coerce"),
            "Ka_Ks": pd.to_numeric(row.get("Ka/Ks"), errors="coerce"),
            "source_file": path.name,
        }
    )

kaks = pd.DataFrame(rows)
pair_to_info = {
    f"{r.parent}-{r.child}": (r.two_method_agreement, r.duplication_type)
    for r in traj.itertuples(index=False)
}
pair_to_info.update({
    f"{r.child}-{r.parent}": (r.two_method_agreement, r.duplication_type)
    for r in traj.itertuples(index=False)
})
kaks["trajectory"] = kaks["pair"].map(lambda pair: pair_to_info.get(pair, (None, None))[0])
kaks["duplication_type"] = kaks["pair"].map(lambda pair: pair_to_info.get(pair, (None, None))[1])
kaks = kaks[kaks["trajectory"].notna() & kaks["Ka_Ks"].notna()].copy()

traj.to_csv(OUT / "two_method_trajectory_tau.tsv", sep="\t", index=False)
kaks.to_csv(OUT / "two_method_trajectory_kaks.tsv", sep="\t", index=False)
print(
    {
        "trajectory_pairs": len(traj),
        "valid_kaks": len(kaks),
        "kaks_by_trajectory": kaks["trajectory"].value_counts().to_dict(),
    }
)




