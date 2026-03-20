#!/usr/bin/env python3

"""
Integrate translocation summaries into editing outcomes. Keeps only the top 2 alleles and translocations above a certain threshold of read counts.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ALLELE_REPAIR_RENAMES = {
    "NHEJ": "NHEJ-Deletion",
    "TemplatedInsertion": "Duplication-Insertion",
    "UnclassifiedInsertion": "Insertion",
    "UnclassifiedInsertion+Substitution": "Insertion+Substitution",
    "TemplatedInsertion+Substitution": "Duplication-Insertion+Substitution",
    "NHEJ+Substitution": "NHEJ-Deletion+Substitution",
    "MMEJ": "MMEJ-Deletion",
    "MMEJ+Substitution": "MMEJ-Deletion+Substitution",
}

GROUPED_REPAIR_RENAMES = {
    "IncorrectHDR": "Imperfect HDR",
    "Duplication-Insertion": "Duplication-insertion",
    "Deletion + Insertion": "Deletion + insertion",
    "PAMdistalNHEJDeletion": "PAM-distal NHEJ deletion",
    "PAMproximalNHEJDeletion": "PAM-proximal NHEJ deletion",
    "BidirectionalNHEJDeletion": "Bidirectional NHEJ deletion",
    "PAMdistalMMEJDeletion": "PAM-distal MMEJ deletion",
    "PAMproximalMMEJDeletion": "PAM-proximal MMEJ deletion",
    "BidirectionalMMEJDeletion": "Bidirectional MMEJ deletion",
    "LargeDeletion": "Intra-chromosomal SV",
    "Translocation": "Inter-chromosomal SV",
}

LEFT_ANCHOR = "CTCCTCAC"  # PAM proximal anchor
RIGHT_ANCHOR = "AGTTGCCATG"  # PAM distal anchor
OPTIONAL_METADATA_COLUMNS = [
    "Cell_Index",
    "Cell_Barcode",
    "Sample",
    "Sample_Tag",
    "Sample_Name",
    "In_Whitelist",
    "In_Unfiltered_Whitelist",
    "Celltype",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--editing-outcomes", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--target-chromosome", default="")
    parser.add_argument("--translocations", default="")
    return parser.parse_args()


def get_max_index(columns: list[str], prefix: str) -> int:
    indices = []
    for column in columns:
        if column.startswith(prefix):
            suffix = column.removeprefix(prefix)
            if suffix.isdigit():
                indices.append(int(suffix))
    return max(indices, default=0)


def add_true_percentages(df: pd.DataFrame) -> pd.DataFrame:
    # Add allele-level percentages relative to the total read and UMI counts for each cell barcode row.
    max_alleles = get_max_index(list(df.columns), "Allele")
    if max_alleles == 0:
        return df

    for allele_index in range(1, max_alleles + 1):
        rc_col = f"rc_{allele_index}"
        umi_col = f"UMIcount_{allele_index}"
        rc_perc_col = f"true_rc_perc_{allele_index}"
        umi_perc_col = f"UMIcount_perc_{allele_index}"

        if rc_col in df.columns:
            df[rc_perc_col] = df[rc_col] / df["total_rc"] * 100
        if umi_col in df.columns:
            df[umi_perc_col] = df[umi_col] / df["total_UMIcount"] * 100

    return df


def top_alleles(df: pd.DataFrame, threshold: float = 10.0, top_n: int = 2) -> pd.DataFrame:
    # Keep only the dominant allele calls per cell to simplify downstream reporting and translocation merging.
    max_alleles = get_max_index(list(df.columns), "Allele")
    metadata_cols = [
        column
        for column in [
            "bc",
            "allele_count",
            "total_rc",
            "total_UMIcount",
            *OPTIONAL_METADATA_COLUMNS,
        ]
        if column in df.columns
    ]
    rows: list[dict] = []

    for _, row in df.iterrows():
        allele_entries = []
        for allele_index in range(1, max_alleles + 1):
            umi_col = f"UMIcount_{allele_index}"
            if umi_col not in df.columns or pd.isna(row.get(umi_col)):
                continue
            entry = {
                "Allele": row.get(f"Allele{allele_index}"),
                "rc": row.get(f"rc_{allele_index}"),
                "rc_perc": row.get(f"rc_perc_{allele_index}"),
                "true_rc_perc": row.get(f"true_rc_perc_{allele_index}"),
                "Sequence": row.get(f"Sequence_{allele_index}"),
                "UMIcount": row.get(umi_col),
                "UMIcount_perc": row.get(f"UMIcount_perc_{allele_index}"),
                "HDist_HDRBC": row.get(f"HDist_HDRBC_{allele_index}"),
                "AlnRefSeq": row.get(f"AlnRefSeq_{allele_index}"),
                "RepairOutcome": row.get(f"RepairOutcome_{allele_index}"),
            }
            if pd.notna(entry["rc_perc"]) and entry["rc_perc"] > threshold:
                allele_entries.append(entry)

        if not allele_entries:
            continue

        top_entries = sorted(allele_entries, key=lambda entry: entry["rc_perc"], reverse=True)[:top_n]
        output_row = {column: row.get(column) for column in metadata_cols}
        output_row["allele_count"] = len(top_entries)

        for allele_index in range(1, top_n + 1):
            if allele_index <= len(top_entries):
                entry = top_entries[allele_index - 1]
                output_row[f"Allele{allele_index}"] = entry["Allele"]
                output_row[f"rc_{allele_index}"] = entry["rc"]
                output_row[f"rc_perc_{allele_index}"] = entry["rc_perc"]
                output_row[f"true_rc_perc_{allele_index}"] = entry["true_rc_perc"]
                output_row[f"Sequence_{allele_index}"] = entry["Sequence"]
                output_row[f"UMIcount_{allele_index}"] = entry["UMIcount"]
                output_row[f"UMIcount_perc_{allele_index}"] = entry["UMIcount_perc"]
                output_row[f"HDist_HDRBC_{allele_index}"] = entry["HDist_HDRBC"]
                output_row[f"AlnRefSeq_{allele_index}"] = entry["AlnRefSeq"]
                output_row[f"RepairOutcome_{allele_index}"] = entry["RepairOutcome"]
            else:
                output_row[f"Allele{allele_index}"] = None
                output_row[f"rc_{allele_index}"] = None
                output_row[f"rc_perc_{allele_index}"] = None
                output_row[f"true_rc_perc_{allele_index}"] = None
                output_row[f"Sequence_{allele_index}"] = None
                output_row[f"UMIcount_{allele_index}"] = None
                output_row[f"UMIcount_perc_{allele_index}"] = None
                output_row[f"HDist_HDRBC_{allele_index}"] = None
                output_row[f"AlnRefSeq_{allele_index}"] = None
                output_row[f"RepairOutcome_{allele_index}"] = None

        rows.append(output_row)

    return pd.DataFrame(rows)


def outcome_frequencies(df: pd.DataFrame, prefix: str, label: str) -> pd.DataFrame:
    outcome_columns = [column for column in df.columns if column.startswith(prefix)]
    if not outcome_columns:
        return pd.DataFrame(columns=["OutcomeType", "Outcome", "Frequency", "Percentage"])

    values = pd.concat([df[column] for column in outcome_columns]).dropna()
    if values.empty:
        return pd.DataFrame(columns=["OutcomeType", "Outcome", "Frequency", "Percentage"])

    frequencies = values.value_counts()
    percentages = frequencies / frequencies.sum() * 100
    return pd.DataFrame(
        {
            "OutcomeType": label,
            "Outcome": frequencies.index,
            "Frequency": frequencies.values,
            "Percentage": percentages.values,
        }
    )


def load_translocations(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    df = pd.read_csv(path).fillna(0)
    if df.empty or "bc" not in df.columns:
        return None
    df["bc"] = df["bc"].astype(str)
    return df


def add_translocation_percentages(df: pd.DataFrame) -> pd.DataFrame:
    max_locations = get_max_index(list(df.columns), "location_")
    for location_index in range(1, max_locations + 1):
        umi_col = f"umi_count_loc_{location_index}"
        rc_col = f"read_count_loc_{location_index}"
        if umi_col in df.columns:
            df[f"translocation_UMIcount_perc_{location_index}"] = df[umi_col] / df["total_UMIcount"] * 100
        if rc_col in df.columns:
            df[f"translocation_rc_perc_{location_index}"] = df[rc_col] / df["total_rc"] * 100
    return df


def top_translocations(df: pd.DataFrame, threshold: float = 10.0, top_n: int = 2) -> pd.DataFrame:
    max_locations = get_max_index(list(df.columns), "location_")
    metadata_cols = [
        column
        for column in [
            "bc",
            "total_rc",
            "total_UMIcount",
            *OPTIONAL_METADATA_COLUMNS,
        ]
        if column in df.columns
    ]
    rows: list[dict] = []

    for _, row in df.iterrows():
        location_entries = []
        for location_index in range(1, max_locations + 1):
            umi_col = f"umi_count_loc_{location_index}"
            if umi_col not in df.columns or pd.isna(row.get(umi_col)):
                continue
            entry = {
                "location": row.get(f"location_{location_index}"),
                "umi_count_loc": row.get(umi_col),
                "translocation_UMIcount_perc": row.get(f"translocation_UMIcount_perc_{location_index}"),
                "read_count_loc": row.get(f"read_count_loc_{location_index}"),
                "translocation_rc_perc": row.get(f"translocation_rc_perc_{location_index}"),
                "read_seq_loc": row.get(f"read_seq_loc_{location_index}"),
            }
            if (
                pd.notna(entry["translocation_rc_perc"])
                and pd.notna(entry["translocation_UMIcount_perc"])
                and entry["translocation_rc_perc"] > threshold
                and entry["translocation_UMIcount_perc"] > threshold
            ):
                location_entries.append(entry)

        if not location_entries:
            continue

        top_entries = sorted(location_entries, key=lambda entry: entry["translocation_rc_perc"], reverse=True)[:top_n]
        output_row = {column: row.get(column) for column in metadata_cols}
        output_row["location_count"] = len(top_entries)

        for location_index in range(1, top_n + 1):
            if location_index <= len(top_entries):
                entry = top_entries[location_index - 1]
                output_row[f"location_{location_index}"] = entry["location"]
                output_row[f"umi_count_loc_{location_index}"] = entry["umi_count_loc"]
                output_row[f"translocation_UMIcount_perc_{location_index}"] = entry["translocation_UMIcount_perc"]
                output_row[f"read_count_loc_{location_index}"] = entry["read_count_loc"]
                output_row[f"translocation_rc_perc_{location_index}"] = entry["translocation_rc_perc"]
                output_row[f"read_seq_loc_{location_index}"] = entry["read_seq_loc"]
            else:
                output_row[f"location_{location_index}"] = None
                output_row[f"umi_count_loc_{location_index}"] = None
                output_row[f"translocation_UMIcount_perc_{location_index}"] = None
                output_row[f"read_count_loc_{location_index}"] = None
                output_row[f"translocation_rc_perc_{location_index}"] = None
                output_row[f"read_seq_loc_{location_index}"] = None

        rows.append(output_row)

    if rows:
        return pd.DataFrame(rows)

    empty_columns = [*metadata_cols, "location_count"]
    for location_index in range(1, top_n + 1):
        empty_columns.extend(
            [
                f"location_{location_index}",
                f"umi_count_loc_{location_index}",
                f"translocation_UMIcount_perc_{location_index}",
                f"read_count_loc_{location_index}",
                f"translocation_rc_perc_{location_index}",
                f"read_seq_loc_{location_index}",
            ]
        )
    return pd.DataFrame(columns=empty_columns)


def merge_alleles_and_translocations(
    allele_df: pd.DataFrame,
    translocation_df: pd.DataFrame | None,
) -> pd.DataFrame:
    # Combine editing outcomes with remapped translocation summaries while keeping total counts consistent across both data sources.
    if translocation_df is None or translocation_df.empty:
        merged = allele_df.copy()
        merged["location_count"] = 0
        return merged

    translocation_totals = translocation_df[["bc", "total_umis", "total_reads"]].copy()
    if not allele_df.empty:
        totals = allele_df[["bc", "total_rc", "total_UMIcount"]].merge(translocation_totals, on="bc", how="outer")
        totals["total_rc"] = totals["total_rc"].fillna(0) + totals["total_reads"].fillna(0)
        totals["total_UMIcount"] = totals["total_UMIcount"].fillna(0) + totals["total_umis"].fillna(0)
        totals = totals[["bc", "total_rc", "total_UMIcount"]]
    else:
        totals = translocation_totals.rename(columns={"total_reads": "total_rc", "total_umis": "total_UMIcount"})

    translocations_with_totals = translocation_df.merge(totals, on="bc", how="left", suffixes=("", "_updated"))
    if "total_rc_updated" in translocations_with_totals.columns:
        translocations_with_totals["total_rc"] = translocations_with_totals["total_rc_updated"]
        translocations_with_totals["total_UMIcount"] = translocations_with_totals["total_UMIcount_updated"]
        translocations_with_totals = translocations_with_totals.drop(
            columns=["total_rc_updated", "total_UMIcount_updated"]
        )
    translocations_with_totals = add_translocation_percentages(translocations_with_totals)
    top_locations = top_translocations(translocations_with_totals)
    if top_locations.empty or "bc" not in top_locations.columns:
        merged = allele_df.copy()
        merged["location_count"] = 0
        return merged

    renamed_top_locations = top_locations.rename(
        columns={
            "total_rc": "translocation_total_rc",
            "total_UMIcount": "translocation_total_UMIcount",
        }
    )
    merged = allele_df.merge(renamed_top_locations, on="bc", how="outer")

    if "translocation_total_rc" in merged.columns:
        merged["total_rc"] = merged["translocation_total_rc"].combine_first(merged.get("total_rc"))
        merged = merged.drop(columns=["translocation_total_rc"])
    if "translocation_total_UMIcount" in merged.columns:
        merged["total_UMIcount"] = merged["translocation_total_UMIcount"].combine_first(merged.get("total_UMIcount"))
        merged = merged.drop(columns=["translocation_total_UMIcount"])

    merged["location_count"] = merged.get("location_count", 0).fillna(0)
    return merged


def normalize_outcomes(df: pd.DataFrame, target_chromosome: str) -> pd.DataFrame:
    df = df.copy()
    for allele_index in range(1, 3):
        repair_col = f"RepairOutcome_{allele_index}"
        if repair_col in df.columns:
            df[repair_col] = df[repair_col].replace(ALLELE_REPAIR_RENAMES)

    def classify_location(location: str | float | None) -> str | None:
        if pd.isna(location):
            return None
        location = str(location)
        if target_chromosome and location.startswith(f"{target_chromosome}:"):
            return "LargeDeletion"
        if location.startswith("chr"):
            return "Translocation"
        return location

    def assign_row(row: pd.Series) -> pd.Series:
        location_pairs = [
            (row.get("location_1"), classify_location(row.get("location_1"))),
            (row.get("location_2"), classify_location(row.get("location_2"))),
        ]
        existing = [
            row.get("RepairOutcome_1") if pd.notna(row.get("RepairOutcome_1")) else None,
            row.get("RepairOutcome_2") if pd.notna(row.get("RepairOutcome_2")) else None,
        ]
        open_slots = [index for index, value in enumerate(existing) if value is None]

        for raw_location, label in location_pairs:
            if label is None or pd.isna(raw_location):
                continue
            if open_slots:
                slot = open_slots.pop(0) + 1
                row[f"RepairOutcome_{slot}"] = label
                row[f"Allele{slot}"] = raw_location
            elif pd.isna(row.get("RepairOutcome_2")):
                row["RepairOutcome_2"] = label
                row["Allele2"] = raw_location
        return row

    return df.apply(assign_row, axis=1)


def assign_deletion_group(sequence: object, outcome: object) -> object:
    if pd.isna(sequence) or pd.isna(outcome) or outcome == "":
        return outcome

    sequence = str(sequence)
    outcome = str(outcome)

    if "NHEJ" in outcome:
        if LEFT_ANCHOR not in sequence and RIGHT_ANCHOR not in sequence:
            return "BidirectionalNHEJDeletion"
        if LEFT_ANCHOR not in sequence:
            return "PAMproximalNHEJDeletion"
        if RIGHT_ANCHOR not in sequence:
            return "PAMdistalNHEJDeletion"
    if "MMEJ" in outcome:
        if LEFT_ANCHOR not in sequence and RIGHT_ANCHOR not in sequence:
            return "BidirectionalMMEJDeletion"
        if LEFT_ANCHOR not in sequence:
            return "PAMproximalMMEJDeletion"
        if RIGHT_ANCHOR not in sequence:
            return "PAMdistalMMEJDeletion"
    return outcome


def group_repair_outcome(outcome: object, sequence: object) -> object:
    if pd.isna(outcome):
        return outcome

    outcome_str = str(outcome)
    if "NHEJ" in outcome_str or "MMEJ" in outcome_str:
        outcome_str = str(assign_deletion_group(sequence, outcome_str))
    grouped = ALLELE_REPAIR_RENAMES.get(outcome_str, outcome_str)
    return GROUPED_REPAIR_RENAMES.get(grouped, grouped)


def classify_repair_outcome(outcome: object) -> object:
    if pd.isna(outcome):
        return outcome

    outcome_lower = str(outcome).lower()
    if "gc" in outcome_lower:
        return "Insertion"
    if "deletion" in outcome_lower and "insertion" in outcome_lower:
        return "Deletion + Insertion"
    if "mmej" in outcome_lower:
        return "MMEJ-Deletion"
    if "nhej" in outcome_lower:
        return "NHEJ-Deletion"
    if "substitution_" in outcome_lower:
        return "Substitution"
    if "duplication-insertion" in outcome_lower:
        return "Insertion"
    if "itrinsertion" in outcome_lower:
        return "ITRInsertion"
    if "repeatedtemplatedinsertion" in outcome_lower:
        return "Insertion"
    if "insertion" in outcome_lower:
        return "Insertion"
    if "incompletehdr" in outcome_lower or "incorrecthdr" in outcome_lower or "imperfect hdr" in outcome_lower:
        return "HDR"
    if "translocation" in outcome_lower or "sv" in outcome_lower or "largedeletion" in outcome_lower:
        return "SV"
    return outcome


def add_grouped_and_class_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Add higher-order grouped labels and broad repair classes next to the raw
    # outcome columns for easier plotting and downstream summarization.
    df = df.copy()
    if "RepairOutcome_1" in df.columns and "RepairOutcome_2" in df.columns:
        df["RepairOutcome_2"] = df["RepairOutcome_2"].fillna(df["RepairOutcome_1"])
    if "Allele1" in df.columns and "Allele2" in df.columns:
        df["Allele2"] = df["Allele2"].fillna(df["Allele1"])
    for prefix in ("Sequence_", "AlnRefSeq_", "HDist_HDRBC_"):
        source_col = f"{prefix}1"
        target_col = f"{prefix}2"
        if source_col in df.columns and target_col in df.columns:
            df[target_col] = df[target_col].fillna(df[source_col])

    for slot in (1, 2):
        repair_col = f"RepairOutcome_{slot}"
        grouped_col = f"GroupedRepairOutcome_{slot}"
        class_col = f"RepairClass_{slot}"
        if repair_col in df.columns:
            df[grouped_col] = df.apply(
                lambda row: group_repair_outcome(row[repair_col], row.get(f"Sequence_{slot}")),
                axis=1,
            )
            df[class_col] = df[grouped_col].apply(classify_repair_outcome)

    ordered_columns = []
    for column in df.columns:
        if column in {"GroupedRepairOutcome_1", "GroupedRepairOutcome_2", "RepairClass_1", "RepairClass_2"}:
            continue
        ordered_columns.append(column)
        if column == "RepairOutcome_1" and "GroupedRepairOutcome_1" in df.columns:
            ordered_columns.append("GroupedRepairOutcome_1")
            ordered_columns.append("RepairClass_1")
        if column == "RepairOutcome_2" and "GroupedRepairOutcome_2" in df.columns:
            ordered_columns.append("GroupedRepairOutcome_2")
            ordered_columns.append("RepairClass_2")

    return df[ordered_columns]


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    allele_df = pd.read_excel(args.editing_outcomes)
    allele_df["bc"] = allele_df["bc"].astype(str)
    allele_df = add_true_percentages(allele_df)
    top_allele_df = top_alleles(allele_df)

    translocation_df = load_translocations(Path(args.translocations)) if args.translocations else None
    merged_df = merge_alleles_and_translocations(top_allele_df, translocation_df)
    merged_df = normalize_outcomes(merged_df, args.target_chromosome)

    merged_path = output_dir / "EditingOutcomesWithTranslocations.csv"
    merged_df.to_csv(merged_path, index=False)

    filtered_df = merged_df.copy()
    for allele_index in range(1, 3):
        if f"UMIcount_{allele_index}" in filtered_df.columns and f"total_UMIcount" in filtered_df.columns:
            filtered_df[f"UMIcount_perc_{allele_index}"] = filtered_df[f"UMIcount_{allele_index}"] / filtered_df["total_UMIcount"] * 100
        if f"rc_{allele_index}" in filtered_df.columns and f"total_rc" in filtered_df.columns:
            filtered_df[f"rc_perc_{allele_index}"] = filtered_df[f"rc_{allele_index}"] / filtered_df["total_rc"] * 100
        if f"umi_count_loc_{allele_index}" in filtered_df.columns and f"total_UMIcount" in filtered_df.columns:
            filtered_df[f"translocation_UMIcount_perc_{allele_index}"] = filtered_df[f"umi_count_loc_{allele_index}"] / filtered_df["total_UMIcount"] * 100
        if f"read_count_loc_{allele_index}" in filtered_df.columns and f"total_rc" in filtered_df.columns:
            filtered_df[f"translocation_rc_perc_{allele_index}"] = filtered_df[f"read_count_loc_{allele_index}"] / filtered_df["total_rc"] * 100

    filtered_df = add_grouped_and_class_columns(filtered_df)
    pd.concat(
        [
            outcome_frequencies(filtered_df, "RepairOutcome_", "RepairOutcome"),
            outcome_frequencies(filtered_df, "GroupedRepairOutcome_", "GroupedRepairOutcome"),
            outcome_frequencies(filtered_df, "RepairClass_", "RepairClass"),
        ],
        ignore_index=True,
    ).to_csv(output_dir / "EditingOutcomeFrequencies.csv", index=False)
    filtered_df.to_csv(output_dir / "FilteredEditingOutcomesWithTranslocations.csv", index=False)


if __name__ == "__main__":
    main()
