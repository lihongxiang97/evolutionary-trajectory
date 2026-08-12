#!/usr/bin/env python3
"""Audit source-data integrity and the expected manuscript-figure coverage."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source_data"
INDEX = SOURCE / "SOURCE_DATA_INDEX.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    expected = [
        *(SOURCE / "main_figures" / name for name in (
            "figure1", "figure2", "figure3", "figure4", "figures5_7", "figure8"
        )),
        SOURCE / "supplementary_figures",
        ROOT / "analysis_results" / "3D_duplicate_mechanisms",
        ROOT / "analysis_results" / "Cross_Haplotype_Robustness",
        ROOT / "analysis_results" / "Risk_Audit",
    ]
    missing = [str(path.relative_to(ROOT)) for path in expected if not path.is_dir()]
    if missing:
        raise SystemExit(f"Missing expected directories: {missing}")

    with INDEX.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    failures: list[str] = []
    for row in rows:
        path = ROOT / row["path"]
        if not path.is_file():
            failures.append(f"missing: {row['path']}")
            continue
        if path.stat().st_size != int(row["size_bytes"]):
            failures.append(f"size: {row['path']}")
        if sha256(path) != row["sha256"]:
            failures.append(f"sha256: {row['path']}")
        if path.suffix == ".json":
            json.loads(path.read_text(encoding="utf-8"))
        elif path.suffix == ".gz":
            with gzip.open(path, "rb") as handle:
                handle.read(1)
        elif path.suffix in {".tsv", ".bed", ".bedgraph", ".txt"}:
            with path.open("rb") as handle:
                handle.read(1)

    if failures:
        raise SystemExit("Source-data audit failed:\n" + "\n".join(failures))
    print(f"PASS: {len(rows)} indexed source-data files; all expected figure/result groups present")


if __name__ == "__main__":
    main()
