#!/usr/bin/env python3
"""Regenerate the source-data inventory with SHA-256 checksums."""
from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source_data"
OUTPUT = SOURCE / "SOURCE_DATA_INDEX.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


rows = []
for path in sorted(p for p in SOURCE.rglob("*") if p.is_file() and p != OUTPUT):
    rows.append((path.relative_to(ROOT).as_posix(), path.stat().st_size, sha256(path)))

with OUTPUT.open("w", encoding="utf-8", newline="\n") as handle:
    handle.write("path\tsize_bytes\tsha256\n")
    for relative, size, digest in rows:
        handle.write(f"{relative}\t{size}\t{digest}\n")

print(f"Indexed {len(rows)} source-data files: {OUTPUT}")
