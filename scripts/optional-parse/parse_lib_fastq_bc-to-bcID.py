#!/usr/bin/env python3

"""Rewrite PARSE FASTQ headers with barcode IDs from local lookup tables."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def hamming_distance(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("Barcode lengths are not equal.")
    return sum(base_left != base_right for base_left, base_right in zip(left, right))


def load_bc_map(barcode_file: Path) -> dict[str, tuple[str, str]]:
    barcode_map: dict[str, tuple[str, str]] = {}
    with barcode_file.open() as handle:
      for line in handle:
            parts = line.strip().split(",")
            barcode_map[parts[1]] = (parts[0], parts[4])
    return barcode_map


def find_closest_bc(barcode: str, barcode_map: dict[str, tuple[str, str]]) -> tuple[str, str]:
    close_barcodes = [candidate for candidate in barcode_map if hamming_distance(barcode, candidate) == 1]
    if len(close_barcodes) == 1:
        return barcode_map[close_barcodes[0]]
    return "-", "-"


def match_ids(header: str, bc_map_24: dict[str, tuple[str, str]], bc_map_v1: dict[str, tuple[str, str]]) -> str:
    parts = header.split("_")
    bc1 = parts[1][0:8]
    bc2 = parts[1][8:16]
    bc3 = parts[1][16:24]

    bc1_id, sample_type_1 = bc_map_24.get(bc1, ("-", "-"))
    bc2_id, sample_type_2 = bc_map_v1.get(bc2, ("-", "-"))
    bc3_id, sample_type_3 = bc_map_v1.get(bc3, ("-", "-"))

    if bc1_id == "-":
        bc1_id, sample_type_1 = find_closest_bc(bc1, bc_map_24)
    if bc2_id == "-":
        bc2_id, sample_type_2 = find_closest_bc(bc2, bc_map_v1)
    if bc3_id == "-":
        bc3_id, sample_type_3 = find_closest_bc(bc3, bc_map_v1)

    return (
        f"{parts[0]}_{bc1_id}:{bc2_id}:{bc3_id}_{bc1}:{bc2}:{bc3}_"
        f"{sample_type_1}_{'_'.join(parts[2:])}"
    )


def process_fastq(
    infile: Path,
    outfile: Path,
    bc_map_24: dict[str, tuple[str, str]],
    bc_map_v1: dict[str, tuple[str, str]],
) -> None:
    with gzip.open(infile, "rt") as fin, gzip.open(outfile, "wt") as fout:
        while True:
            lines = [next(fin, "").strip() for _ in range(4)]
            if not lines[0]:
                break
            lines[0] = match_ids(lines[0], bc_map_24, bc_map_v1)
            fout.write("\n".join(lines) + "\n")


def main() -> None:
    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2])
    parsekit = sys.argv[3] if len(sys.argv) > 3 else ""

    bcfile24 = SCRIPT_DIR / ("bc_data_n24_v4.csv" if parsekit == "mini" else "bc_data_n198_v5.csv")
    bcfilev1 = SCRIPT_DIR / "bc_data_v1.csv"

    process_fastq(infile, outfile, load_bc_map(bcfile24), load_bc_map(bcfilev1))


if __name__ == "__main__":
    main()
