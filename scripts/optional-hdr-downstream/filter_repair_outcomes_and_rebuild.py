#!/usr/bin/env python3

"""
Filter UMIs from an existing RepairOutcomeTable.csv and rebuild downstream
editing-outcome summary files without rerunning CRISPResso or remapping.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "required-core"))

from bdr_utils import load_bdr_sample_tags, normalize_bdr_barcode_series
from dro_target_config import load_dro_config
from integrate_translocations_and_filter_outcomes import (
    add_grouped_and_class_columns,
    add_true_percentages,
    load_translocations,
    merge_alleles_and_translocations,
    normalize_outcomes,
    outcome_frequencies,
    top_alleles,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repair-table", required=True)
    parser.add_argument("--outcomes-list", required=True)
    parser.add_argument("--cbc-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--libtype", default="10X")
    parser.add_argument("--target-chromosome", default="")
    parser.add_argument("--translocations", default="")
    return parser.parse_args()


def hamming_distance(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def nonempty_string(value: object) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip() not in {"", "0", "0.0", "nan", "None"}


def dominant_mod(mods: object) -> str | None:
    if not nonempty_string(mods):
        return None
    first_entry = str(mods).split(";")[0]
    if not first_entry:
        return None
    return first_entry.split("_", 1)[0]


def load_outcome_filter(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def parse_mod_entry(entry: str) -> tuple[str, int, str, int, str] | None:
    parts = entry.split("_", 4)
    if len(parts) != 5:
        return None
    mod, count, modseq, modseq_count, refalnseq = parts
    try:
        return mod, int(count), modseq, int(modseq_count), refalnseq
    except ValueError:
        return None


def parse_hdr_entry(entry: str) -> tuple[str, int, str, int, int] | None:
    parts = entry.split("_", 4)
    if len(parts) != 5:
        return None
    hdr_barcode, count, hdrbc_seq, hdrbc_seq_count, hdrbc_hd = parts
    try:
        return hdr_barcode, int(count), hdrbc_seq, int(hdrbc_seq_count), int(hdrbc_hd)
    except ValueError:
        return None


def parse_wt_entry(entry: object) -> tuple[str, int] | None:
    if not nonempty_string(entry):
        return None
    parts = str(entry).split("_", 1)
    if len(parts) != 2:
        return None
    seq, count = parts
    try:
        return seq, int(count)
    except ValueError:
        return None


def find_first_gap_in_alnseq(aligned_seq: str) -> int | None:
    for index, char in enumerate(aligned_seq):
        if char == "-":
            return index
    return None


def extract_deletion_surrounding_bases(aligned_seq: str, gap_pos: int, mod_length: int) -> tuple[str, str]:
    start_pos = max(0, gap_pos - 25)
    left_bases = aligned_seq[start_pos:gap_pos].replace("-", "")
    gap_end_pos = gap_pos + mod_length
    right_bases = aligned_seq[gap_end_pos:gap_end_pos + 25].replace("-", "")
    return left_bases, right_bases


def extract_deleted_bases(aligned_refseq: str, aligned_seq: str, mod_length: int) -> str | None:
    gap_pos = find_first_gap_in_alnseq(aligned_seq)
    if gap_pos is None:
        return None
    return aligned_refseq[gap_pos:gap_pos + mod_length]


def find_longest_microhomology(
    deleted_seq: str,
    left_part: str,
    right_part: str,
    min_homology_length: int = 2,
    tolerance: int = 2,
) -> tuple[str, str, str] | None:
    def find_longest_start(seq1: str, seq2: str) -> str | None:
        longest = ""
        for length in range(min_homology_length, min(len(seq1), len(seq2)) + 1):
            for start1 in range(min(tolerance + 1, len(seq1) - length + 1)):
                for start2 in range(min(tolerance + 1, len(seq2) - length + 1)):
                    if seq1[start1:start1 + length] == seq2[start2:start2 + length] and length > len(longest):
                        longest = seq1[start1:start1 + length]
        return longest or None

    def find_longest_end(seq1: str, seq2: str) -> str | None:
        longest = ""
        for length in range(min_homology_length, min(len(seq1), len(seq2)) + 1):
            for end1 in range(len(seq1) - tolerance, len(seq1) - length + 1):
                for end2 in range(len(seq2) - tolerance, len(seq2) - length + 1):
                    if seq1[end1:end1 + length] == seq2[end2:end2 + length] and length > len(longest):
                        longest = seq1[end1:end1 + length]
        return longest or None

    candidates = [
        (find_longest_start(deleted_seq, left_part), "start-to-start", "left"),
        (find_longest_start(deleted_seq, right_part), "start-to-start", "right"),
        (find_longest_end(deleted_seq, left_part), "end-to-end", "left"),
        (find_longest_end(deleted_seq, right_part), "end-to-end", "right"),
    ]

    longest: tuple[str, str, str] | None = None
    for candidate, mh_type, source in candidates:
        if candidate and (longest is None or len(candidate) > len(longest[0])):
            longest = (candidate, mh_type, source)
    return longest


def has_near_match(inserted_bases: str, reference_sequence: str, max_mismatches: int = 0) -> bool:
    for index in range(len(reference_sequence) - len(inserted_bases) + 1):
        reference = reference_sequence[index:index + len(inserted_bases)]
        if hamming_distance(inserted_bases, reference) <= max_mismatches:
            return True
    return False


def find_gap_in_refseq(aligned_refseq: str) -> int | None:
    gap_positions = [index for index, char in enumerate(aligned_refseq) if char == "-"]
    return gap_positions[0] if gap_positions else None


def extract_surrounding_bases(aligned_refseq: str, aligned_seq: str, mod_length: int) -> tuple[str, str] | None:
    gap_pos = find_gap_in_refseq(aligned_refseq)
    aligned_seq = aligned_seq.rstrip("-")
    if gap_pos is None:
        return None
    start_pos = max(0, gap_pos - mod_length)
    end_pos = min(len(aligned_refseq), gap_pos + mod_length + mod_length)
    surrounding_ref_bases = aligned_refseq[start_pos:end_pos].replace("-", "")
    inserted_bases = str(aligned_seq[gap_pos:gap_pos + mod_length])
    return surrounding_ref_bases, inserted_bases


def find_repeating_base_element(sequence: str) -> str:
    for length in range(1, len(sequence) // 2 + 1):
        if len(sequence) % length == 0:
            base_element = sequence[:length]
            if base_element * (len(sequence) // length) == sequence:
                return base_element
    return sequence


def classify_insertion(
    mod: str,
    aligned_refseq: str,
    aligned_seq: str,
    insertion_type: str,
    cins_hdrbc: str | None,
) -> str | bool:
    aav_itr = (
        "AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAG"
    )
    aav_itr_rev = (
        "CTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT"
    )
    gc_set = {"G", "C"}

    mod_length = int(mod[1:])
    extracted = extract_surrounding_bases(aligned_refseq, aligned_seq, mod_length)
    if extracted is None:
        return False
    surrounding_ref_bases, inserted_bases = extracted

    if insertion_type == "cINS" and cins_hdrbc:
        if len(inserted_bases) >= 3 and inserted_bases in cins_hdrbc:
            return "HDR" if len(inserted_bases) == 12 else f"IncompleteHDR-{len(inserted_bases)}"

    if inserted_bases in surrounding_ref_bases:
        return "TemplatedInsertion"
    if len(inserted_bases) > 1 and all(base in gc_set for base in inserted_bases):
        return "GCInsertion"
    if len(inserted_bases) >= 8:
        max_mismatches = 0 if len(inserted_bases) <= 20 else (3 if len(inserted_bases) <= 40 else 4)
        if has_near_match(inserted_bases, aav_itr, max_mismatches) or has_near_match(inserted_bases, aav_itr_rev, max_mismatches):
            return "ITRInsertion"
    base_element = find_repeating_base_element(inserted_bases)
    if base_element in surrounding_ref_bases:
        return "RepeatedTemplatedInsertion"
    return False


def classify_repair_outcome(
    allele: object,
    sequence: object,
    aln_ref_seq: object,
    hdist_hdrbc: object,
    hdrbc_len: int,
    hdr_anchor: str,
    insertion_type: str,
    cins_hdrbc: str | None,
) -> str:
    if pd.isna(allele) or allele in {"", None}:
        return ""

    allele = str(allele)
    cleaned_sequence = str(sequence).replace("-", "")
    hamming_distance_value = int(float(hdist_hdrbc)) if pd.notna(hdist_hdrbc) and str(hdist_hdrbc) != "" else 0

    if len(allele) == hdrbc_len:
        return "HDR" if hamming_distance_value >= (hdrbc_len - 1) else f"IncompleteHDR-{hamming_distance_value}"
    if allele == "WT":
        return "WT"

    match = re.search(r"(D(\d+))?(I(\d+))?(S(\d+))?", allele)
    if not match:
        return ""

    deletion_length = int(match.group(2)) if match.group(2) else 0
    insertion_length = int(match.group(4)) if match.group(4) else 0
    substitution_length = int(match.group(6)) if match.group(6) else 0

    if insertion_length > 0 and deletion_length > 0 and substitution_length > 0:
        return "IncorrectHDR" if hdr_anchor in cleaned_sequence else "Deletion+Insertion+Substitution"
    if insertion_length > 0 and deletion_length == 0 and substitution_length > 0:
        if hdr_anchor in cleaned_sequence:
            return "IncorrectHDR"
        insertion_class = classify_insertion(f"I{insertion_length}", str(aln_ref_seq), str(sequence), insertion_type, cins_hdrbc)
        return f"{insertion_class}+Substitution" if insertion_class else "UnclassifiedInsertion+Substitution"
    if insertion_length > 0 and deletion_length > 0 and substitution_length == 0:
        return "Deletion+Insertion"
    if insertion_length == 0 and deletion_length > 0 and substitution_length == 0:
        deleted_seq = extract_deleted_bases(str(aln_ref_seq), str(sequence), deletion_length)
        if deleted_seq is None:
            return "NHEJ"
        left_part, right_part = extract_deletion_surrounding_bases(str(sequence), find_first_gap_in_alnseq(str(sequence)) or 0, deletion_length)
        result = find_longest_microhomology(deleted_seq, left_part, right_part)
        return "MMEJ" if result and len(result[0]) >= 2 else "NHEJ"
    if insertion_length == 0 and deletion_length > 0 and substitution_length > 0:
        if hdr_anchor in cleaned_sequence:
            return "IncorrectHDR"
        deleted_seq = extract_deleted_bases(str(aln_ref_seq), str(sequence), deletion_length)
        if deleted_seq is None:
            return "UnclassifiedDeletion+Substitution"
        left_part, right_part = extract_deletion_surrounding_bases(str(sequence), find_first_gap_in_alnseq(str(sequence)) or 0, deletion_length)
        result = find_longest_microhomology(deleted_seq, left_part, right_part)
        return "MMEJ+Substitution" if result and len(result[0]) >= 2 else "NHEJ+Substitution"
    if insertion_length == 0 and deletion_length == 0 and substitution_length > 0:
        return "IncorrectHDR" if hdr_anchor in cleaned_sequence else f"Substitution_{substitution_length}"
    if insertion_length > 0 and deletion_length == 0 and substitution_length == 0:
        insertion_class = classify_insertion(allele, str(aln_ref_seq), str(sequence), insertion_type, cins_hdrbc)
        return str(insertion_class) if insertion_class else "UnclassifiedInsertion"
    return ""


def categorize_sized_outcome(allele: object, hdrbc_len: int) -> object:
    if pd.isna(allele):
        return allele
    allele = str(allele)
    if len(allele) == hdrbc_len:
        return "HDR"
    if allele == "WT":
        return "WT"
    if allele.startswith("D"):
        match = re.match(r"D(\d+)", allele)
        if not match:
            return allele
        size = int(match.group(1))
        if size <= 8:
            return "D1-D8"
        if size <= 30:
            return "D9-D30"
        return "D31+"
    if allele.startswith("I"):
        match = re.match(r"I(\d+)", allele)
        if not match:
            return allele
        size = int(match.group(1))
        if size <= 3:
            return "I1-I3"
        if size <= 20:
            return "I4-I20"
        return "I21+"
    return allele


def summarize_hdr_umIs(rotable: pd.DataFrame) -> tuple[dict, dict]:
    umi_hdrbc: dict[tuple[str, str], list[list[object]]] = {}
    cell_hdrbc: dict[str, dict[str, dict[str, object]]] = {}

    for _, row in rotable.iterrows():
        if row["UMI_RepairOutcome"] != "HDR" or not nonempty_string(row.get("hdr_bc")):
            continue
        hdr_read_count = int(row["HDR"])
        true_barcodes = []
        for entry in str(row["hdr_bc"]).split(";"):
            parsed = parse_hdr_entry(entry)
            if parsed is None:
                continue
            hdr_barcode, count, hdrbc_seq, hdrbc_seq_count, hdrbc_hd = parsed
            if hdr_read_count > 0 and count / hdr_read_count > 0.5 and count >= 3:
                true_barcodes.append([hdr_barcode, count, hdr_read_count, hdrbc_seq, hdrbc_seq_count, hdrbc_hd])
        umi_hdrbc[(str(row["bc"]), str(row["umi"]))] = true_barcodes

    for (bc, _umi), values in umi_hdrbc.items():
        if not values:
            continue
        hdrbc = str(values[0][0])
        hdrbc_seq = str(values[0][3])
        hdrbc_hd = int(values[0][5])
        hdrrc = int(values[0][2])
        if bc not in cell_hdrbc:
            cell_hdrbc[bc] = {}
        if hdrbc in cell_hdrbc[bc]:
            cell_hdrbc[bc][hdrbc]["umi_count"] += 1
            cell_hdrbc[bc][hdrbc]["rc"] += hdrrc
        else:
            cell_hdrbc[bc][hdrbc] = {"umi_count": 1, "rc": hdrrc, "seq": hdrbc_seq, "hd": hdrbc_hd}

    for bc, hdrbcs in list(cell_hdrbc.items()):
        if len(hdrbcs) == 2:
            keys = list(hdrbcs.keys())
            if hamming_distance(keys[0], keys[1]) == 1:
                max_key = max(keys, key=lambda key: hdrbcs[key]["umi_count"])
                min_key = min(keys, key=lambda key: hdrbcs[key]["umi_count"])
                hdrbcs[max_key]["umi_count"] += hdrbcs[min_key]["umi_count"]
                hdrbcs[max_key]["rc"] += hdrbcs[min_key]["rc"]
                del hdrbcs[min_key]

    return umi_hdrbc, cell_hdrbc


def summarize_mod_umis(rotable: pd.DataFrame) -> tuple[dict, dict]:
    umi_mod: dict[tuple[str, str], list[list[object]]] = {}
    cell_mod: dict[str, dict[str, dict[str, object]]] = {}

    for _, row in rotable.iterrows():
        if row["UMI_RepairOutcome"] != "NHEJ" or not nonempty_string(row.get("mods")):
            continue
        mods_read_count = max(int(row["Reference_MODIFIED"]), int(row["HDR"]))
        true_mods = []
        for entry in str(row["mods"]).split(";"):
            parsed = parse_mod_entry(entry)
            if parsed is None:
                continue
            mod, count, modseq, modseq_count, refalnseq = parsed
            if mods_read_count > 0 and count / mods_read_count > 0.5 and count >= 3:
                true_mods.append([mod, count, mods_read_count, modseq, modseq_count, refalnseq])
        umi_mod[(str(row["bc"]), str(row["umi"]))] = true_mods

    for (bc, _umi), values in umi_mod.items():
        if not values:
            continue
        mod_entries = values
        if len(mod_entries) > 1:
            modseqs = {str(entry[3]) for entry in mod_entries}
            alnrefseqs = {str(entry[5]) for entry in mod_entries}
            combined_mod = "".join(sorted(str(entry[0]) for entry in mod_entries))
            if len(modseqs) == 1:
                mod_seq = list(modseqs)[0]
                mod_alnrefseq = list(alnrefseqs)[0]
                mod_rc = int(mod_entries[0][2])
                cell_mod.setdefault(bc, {})
                if combined_mod in cell_mod[bc]:
                    cell_mod[bc][combined_mod]["umi_count"] += 1
                    cell_mod[bc][combined_mod]["rc"] += mod_rc
                else:
                    cell_mod[bc][combined_mod] = {
                        "umi_count": 1,
                        "rc": mod_rc,
                        "seq": mod_seq,
                        "alnrefseq": mod_alnrefseq,
                    }
                continue
        for entry in mod_entries:
            mod = str(entry[0])
            mod_seq = str(entry[3])
            mod_alnrefseq = str(entry[5])
            mod_rc = int(entry[2])
            cell_mod.setdefault(bc, {})
            if mod in cell_mod[bc]:
                cell_mod[bc][mod]["umi_count"] += 1
                cell_mod[bc][mod]["rc"] += mod_rc
            else:
                cell_mod[bc][mod] = {
                    "umi_count": 1,
                    "rc": mod_rc,
                    "seq": mod_seq,
                    "alnrefseq": mod_alnrefseq,
                }

    return umi_mod, cell_mod


def summarize_wt_umis(rotable: pd.DataFrame) -> dict[str, dict[str, object]]:
    cell_wt: dict[str, dict[str, object]] = {}
    for _, row in rotable.iterrows():
        if row["UMI_RepairOutcome"] != "WT":
            continue
        parsed = parse_wt_entry(row.get("wt_seq"))
        if parsed is None:
            continue
        wt_seq, _wt_seq_count = parsed
        wt_rc = int(row["Reference_UNMODIFIED"])
        if wt_rc < 3:
            continue
        bc = str(row["bc"])
        if bc in cell_wt:
            cell_wt[bc]["umi_count"] += 1
            cell_wt[bc]["rc"] += wt_rc
        else:
            cell_wt[bc] = {"umi_count": 1, "rc": wt_rc, "seq": wt_seq}
    return cell_wt


def build_cell_alleles(
    rotable: pd.DataFrame,
    cell_hdrbc: dict[str, dict[str, dict[str, object]]],
    cell_mod: dict[str, dict[str, dict[str, object]]],
    cell_wt: dict[str, dict[str, object]],
) -> tuple[dict[str, dict[str, dict[str, object]]], dict[str, int], dict[str, int]]:
    rotable = rotable.copy()
    rotable["totalrc"] = rotable[["HDR", "Reference_UNMODIFIED", "Reference_MODIFIED"]].sum(axis=1)
    cell_totalrc = rotable.groupby("bc")["totalrc"].sum().astype(int).to_dict()
    cell_totalumi = rotable.groupby("bc").size().astype(int).to_dict()

    combined: dict[str, dict[str, dict[str, object]]] = {}
    for bc in sorted(set(rotable["bc"].astype(str))):
        alleles: dict[str, dict[str, object]] = {}
        if bc in cell_hdrbc:
            for allele, values in cell_hdrbc[bc].items():
                rc_perc = values["rc"] / cell_totalrc[bc] * 100 if cell_totalrc[bc] else 0
                alleles[allele] = {**values, "rc_perc": rc_perc}
        if bc in cell_mod:
            for allele, values in cell_mod[bc].items():
                rc_perc = values["rc"] / cell_totalrc[bc] * 100 if cell_totalrc[bc] else 0
                alleles[allele] = {**values, "rc_perc": rc_perc}
        if bc in cell_wt:
            values = cell_wt[bc]
            rc_perc = values["rc"] / cell_totalrc[bc] * 100 if cell_totalrc[bc] else 0
            alleles["WT"] = {**values, "rc_perc": rc_perc}
        if alleles:
            combined[bc] = alleles
    return combined, cell_totalrc, cell_totalumi


def get_metadata_value(metadata_df: pd.DataFrame, barcode: str, *column_names: str, default: str = "Unknown") -> str:
    if barcode not in metadata_df.index:
        return default
    for column_name in column_names:
        if column_name in metadata_df.columns:
            value = metadata_df.loc[barcode, column_name]
            if pd.notna(value):
                return str(value)
    return default


def build_editing_outcomes_dataframe(
    rotable: pd.DataFrame,
    cbc_dir: Path,
    libtype: str,
    hdrbc_len: int,
    hdr_anchor: str,
    insertion_type: str,
    cins_hdrbc: str | None,
) -> pd.DataFrame:
    _, cell_hdrbc = summarize_hdr_umIs(rotable)
    _, cell_mod = summarize_mod_umis(rotable)
    cell_wt = summarize_wt_umis(rotable)
    cell_alleles, cell_totalrc, cell_totalumi = build_cell_alleles(rotable, cell_hdrbc, cell_mod, cell_wt)

    whitelist_df = pd.read_csv(cbc_dir / "cbcs_filtered.csv", low_memory=False)
    whitelist_df["bc"] = whitelist_df["bc"].astype(str)
    whitelist = set(whitelist_df["bc"].tolist())

    whitelist_unfiltered_df = pd.read_csv(cbc_dir / "cbcs_unfiltered.csv", low_memory=False)
    whitelist_unfiltered_df["bc"] = whitelist_unfiltered_df["bc"].astype(str)
    whitelist_unfiltered = set(whitelist_unfiltered_df["bc"].tolist())

    cell_metadata = pd.read_csv(cbc_dir / "cell_metadata.csv", low_memory=False)
    cell_metadata["bc"] = cell_metadata["bc"].astype(str)
    cell_metadata = cell_metadata.set_index("bc")

    sample_tag_metadata = pd.DataFrame()
    if libtype == "BDR":
        sample_tag_metadata = load_bdr_sample_tags()
        if not sample_tag_metadata.empty:
            sample_tag_metadata["bc"] = sample_tag_metadata["bc"].astype(str)
            sample_tag_metadata = sample_tag_metadata.set_index("bc")

    raw_barcodes_by_cell_index: dict[str, str] = {}
    if libtype == "BDR" and "raw_bc" in rotable.columns:
        raw_barcodes_by_cell_index = (
            rotable.loc[rotable["raw_bc"].notna() & (rotable["raw_bc"].astype(str) != ""), ["bc", "raw_bc"]]
            .drop_duplicates(subset=["bc"])
            .assign(bc=lambda df: df["bc"].astype(str), raw_bc=lambda df: df["raw_bc"].astype(str))
            .set_index("bc")["raw_bc"]
            .to_dict()
        )

    rows: list[dict[str, object]] = []
    max_alleles = max((len(alleles) for alleles in cell_alleles.values()), default=0)
    for bc, alleles in cell_alleles.items():
        sorted_alleles = sorted(alleles.items(), key=lambda item: (-int(item[1]["umi_count"]), -int(item[1]["rc"]), item[0]))
        row: dict[str, object] = {
            "bc": bc,
            "allele_count": len(sorted_alleles),
            "total_rc": cell_totalrc.get(bc, 0),
            "total_UMIcount": cell_totalumi.get(bc, 0),
            "Sample": get_metadata_value(cell_metadata, bc, "sample", "Sample"),
            "In_Whitelist": bc in whitelist,
            "In_Unfiltered_Whitelist": bc in whitelist_unfiltered,
            "Celltype": get_metadata_value(cell_metadata, bc, "celltype", "Celltype"),
        }
        if libtype == "BDR":
            row["Cell_Index"] = bc
            row["Cell_Barcode"] = raw_barcodes_by_cell_index.get(bc, "Unknown")
        if not sample_tag_metadata.empty:
            row["Sample_Tag"] = get_metadata_value(sample_tag_metadata, bc, "Sample_Tag")
            row["Sample_Name"] = get_metadata_value(sample_tag_metadata, bc, "Sample_Name")

        detailed_parts = []
        for allele_index, (allele_name, allele_data) in enumerate(sorted_alleles, start=1):
            row[f"Allele{allele_index}"] = allele_name
            row[f"rc_{allele_index}"] = allele_data["rc"]
            row[f"rc_perc_{allele_index}"] = allele_data["rc_perc"]
            row[f"UMIcount_{allele_index}"] = allele_data["umi_count"]
            row[f"Sequence_{allele_index}"] = allele_data.get("seq", "")
            row[f"HDist_HDRBC_{allele_index}"] = allele_data.get("hd", "")
            row[f"AlnRefSeq_{allele_index}"] = allele_data.get("alnrefseq", "")
            row[f"RepairOutcome_{allele_index}"] = classify_repair_outcome(
                allele_name,
                allele_data.get("seq", ""),
                allele_data.get("alnrefseq", ""),
                allele_data.get("hd", ""),
                hdrbc_len,
                hdr_anchor,
                insertion_type,
                cins_hdrbc,
            )
            detailed_parts.append(allele_name)

        for allele_index in range(len(sorted_alleles) + 1, max_alleles + 1):
            row[f"Allele{allele_index}"] = None
            row[f"rc_{allele_index}"] = None
            row[f"rc_perc_{allele_index}"] = None
            row[f"UMIcount_{allele_index}"] = None
            row[f"Sequence_{allele_index}"] = None
            row[f"HDist_HDRBC_{allele_index}"] = None
            row[f"AlnRefSeq_{allele_index}"] = None
            row[f"RepairOutcome_{allele_index}"] = None

        allele_count = len(sorted_alleles)
        total_umicount = int(row["total_UMIcount"])
        if allele_count == 1:
            if total_umicount == 1:
                row["Ploidy"] = "haploid"
                row["DetailedRepairOutcome"] = f"{detailed_parts[0]}_NA"
            else:
                row["Ploidy"] = "diploid"
                row["DetailedRepairOutcome"] = f"{detailed_parts[0]}_{detailed_parts[0]}"
        elif allele_count == 2:
            row["Ploidy"] = "diploid"
            row["DetailedRepairOutcome"] = "_".join(detailed_parts[:2])
        elif allele_count == 3:
            row["Ploidy"] = "triploid"
            row["DetailedRepairOutcome"] = "_".join(detailed_parts[:3])
        else:
            row["Ploidy"] = "polyploid"
            row["DetailedRepairOutcome"] = "_".join(detailed_parts)

        rows.append(row)

    return pd.DataFrame(rows)


def write_frequency_tables(df: pd.DataFrame, output_dir: Path, hdrbc_len: int) -> None:
    if df.empty:
        pd.DataFrame(columns=["Allele", "Frequency"]).to_excel(output_dir / "AlleleFrequencies.xlsx", index=False)
        pd.DataFrame(columns=["Allele", "Frequency"]).to_excel(output_dir / "AlleleFrequencies_sized.xlsx", index=False)
        return

    max_alleles = max(int(col.replace("Allele", "")) for col in df.columns if col.startswith("Allele") and col.replace("Allele", "").isdigit())

    allele_freq = pd.Series(dtype=float)
    for allele_index in range(1, max_alleles + 1):
        column = f"RepairOutcome_{allele_index}"
        if column in df.columns:
            allele_freq = allele_freq.add(df[column].value_counts(), fill_value=0)

    pd.DataFrame({"Allele": allele_freq.index, "Frequency": allele_freq.values}).to_excel(
        output_dir / "AlleleFrequencies.xlsx",
        index=False,
    )

    sized_columns = {}
    for allele_index in range(1, max_alleles + 1):
        allele_column = f"Allele{allele_index}"
        sized_columns[f"Allele{allele_index}_Outcome_sized"] = df[allele_column].apply(lambda value: categorize_sized_outcome(value, hdrbc_len))

    df_sized = pd.concat([df.copy(), pd.DataFrame(sized_columns)], axis=1)

    allele_freq_sized = pd.Series(dtype=float)
    for allele_index in range(1, max_alleles + 1):
        column = f"Allele{allele_index}_Outcome_sized"
        allele_freq_sized = allele_freq_sized.add(df_sized[column].value_counts(), fill_value=0)

    pd.DataFrame({"Allele": allele_freq_sized.index, "Frequency": allele_freq_sized.values}).to_excel(
        output_dir / "AlleleFrequencies_sized.xlsx",
        index=False,
    )
    df_sized.to_excel(output_dir / "EditingOutcomesByCellBarcode.xlsx", index=False, engine="openpyxl")


def main() -> None:
    args = parse_args()
    repair_table_path = Path(args.repair_table).resolve()
    outcomes_path = Path(args.outcomes_list).resolve()
    cbc_dir = Path(args.cbc_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    dro_config = load_dro_config(require_insertion_type=True)
    hdrbc_len = int(dro_config["hdrbc_len"])
    hdr_anchor = dro_config["hdr_anchor"] or "__SCOUT_NO_HDR_ANCHOR__"
    insertion_type = dro_config["insertion_type"]
    cins_hdrbc = dro_config.get("cins_hdrbc")

    outcomes_to_remove = load_outcome_filter(outcomes_path)
    rotable = pd.read_csv(repair_table_path, low_memory=False)
    rotable["bc"] = rotable["bc"].astype(str)
    if args.libtype == "BDR":
        rotable["bc"] = normalize_bdr_barcode_series(rotable["bc"])

    dominant_series = rotable["mods"].apply(dominant_mod) if "mods" in rotable.columns else pd.Series([None] * len(rotable))
    remove_mask = dominant_series.isin(outcomes_to_remove)
    filtered_rotable = rotable.loc[~remove_mask].copy()

    filtered_rotable.to_csv(output_dir / "RepairOutcomeTable.csv", index=False)
    filtered_rotable.groupby("bc")["UMI_RepairOutcome"].value_counts().to_csv(output_dir / "UMIRepairOutcomes.csv")

    editing_df = build_editing_outcomes_dataframe(
        filtered_rotable,
        cbc_dir,
        args.libtype,
        hdrbc_len,
        hdr_anchor,
        insertion_type,
        cins_hdrbc,
    )

    write_frequency_tables(editing_df, output_dir, hdrbc_len)

    allele_df = add_true_percentages(editing_df.copy())
    top_allele_df = top_alleles(allele_df)
    translocation_df = load_translocations(Path(args.translocations)) if args.translocations else None
    merged_df = merge_alleles_and_translocations(top_allele_df, translocation_df)
    merged_df = normalize_outcomes(merged_df, args.target_chromosome)
    merged_df.to_csv(output_dir / "EditingOutcomesWithTranslocations.csv", index=False)

    filtered_df = merged_df.copy()
    for allele_index in range(1, 3):
        if f"UMIcount_{allele_index}" in filtered_df.columns and "total_UMIcount" in filtered_df.columns:
            filtered_df[f"UMIcount_perc_{allele_index}"] = filtered_df[f"UMIcount_{allele_index}"] / filtered_df["total_UMIcount"] * 100
        if f"rc_{allele_index}" in filtered_df.columns and "total_rc" in filtered_df.columns:
            filtered_df[f"rc_perc_{allele_index}"] = filtered_df[f"rc_{allele_index}"] / filtered_df["total_rc"] * 100
        if f"umi_count_loc_{allele_index}" in filtered_df.columns and "total_UMIcount" in filtered_df.columns:
            filtered_df[f"translocation_UMIcount_perc_{allele_index}"] = filtered_df[f"umi_count_loc_{allele_index}"] / filtered_df["total_UMIcount"] * 100
        if f"read_count_loc_{allele_index}" in filtered_df.columns and "total_rc" in filtered_df.columns:
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

    summary = pd.DataFrame(
        [
            {
                "input_rows": len(rotable),
                "removed_rows": int(remove_mask.sum()),
                "kept_rows": len(filtered_rotable),
                "outcomes_filtered": ",".join(sorted(outcomes_to_remove)),
            }
        ]
    )
    summary.to_csv(output_dir / "RepairOutcomeFilterSummary.csv", index=False)


if __name__ == "__main__":
    main()
