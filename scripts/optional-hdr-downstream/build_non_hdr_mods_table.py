#!/usr/bin/env python3

"""
Create an Hdr_mods.csv-compatible table from bc_umi_mod_seq.csv while dropping
all HDR barcode content.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mods-table", required=True, help="Path to bc_umi_mod_seq.csv")
    parser.add_argument("--output", required=True, help="Path to write Hdr_mods.csv")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mods_path = Path(args.mods_table).resolve()
    output_path = Path(args.output).resolve()

    df = pd.read_csv(mods_path, low_memory=False)

    if "HDR" not in df.columns:
        df["HDR"] = ""
    else:
        df["HDR"] = ""

    if "hdr_bc" not in df.columns:
        df["hdr_bc"] = ""
    else:
        df["hdr_bc"] = ""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
