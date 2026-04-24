"""Create a filtered CRISPResso FASTQ that excludes reads matching an anchor-delimited insertion pattern."""

import gzip
import sys


def hdist(a: str, b: str) -> int:
    return len(list(filter(lambda x: ord(x[0]) ^ ord(x[1]), zip(a, b))))


def match_with_mismatch(anchor: str, seq: str, max_mismatches: int = 1):
    for i in range(len(seq) - len(anchor) + 1):
        subseq = seq[i:i + len(anchor)]
        if hdist(anchor, subseq) <= max_mismatches:
            return True, i
    return False, -1


def anchored_gap_length(seq: str, left_anchor: str, right_anchor: str, max_mismatches: int) -> int:
    left_match, left_index = match_with_mismatch(left_anchor, seq, max_mismatches=max_mismatches)
    right_match, right_index = match_with_mismatch(right_anchor, seq, max_mismatches=max_mismatches)
    if left_match and right_match and right_index > left_index:
        return right_index - (left_index + len(left_anchor))
    return -1


def main() -> int:
    if len(sys.argv) != 7:
        print(
            f"Usage: {sys.argv[0]} <input_fastq.gz> <output_fastq.gz> <left_anchor> <right_anchor> <bp_between_anchors> <max_mismatches>",
            file=sys.stderr,
        )
        return 1

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    left_anchor = sys.argv[3]
    right_anchor = sys.argv[4]
    bp_between_anchors = int(sys.argv[5])
    max_mismatches = int(sys.argv[6])
    kept = 0
    removed = 0

    with gzip.open(input_path, "rt") as src, gzip.open(output_path, "wt") as dst:
        while True:
            header = src.readline()
            if not header:
                break
            sequence = src.readline()
            plus = src.readline()
            quality = src.readline()

            if not sequence or not plus or not quality:
                print("Malformed FASTQ: truncated record encountered.", file=sys.stderr)
                return 1

            gap_len = anchored_gap_length(sequence.strip(), left_anchor, right_anchor, max_mismatches)
            if gap_len == bp_between_anchors:
                removed += 1
                continue

            dst.write(header)
            dst.write(sequence)
            dst.write(plus)
            dst.write(quality)
            kept += 1

    print(f"Filtered CRISPResso FASTQ written to {output_path}")
    print(
        f"Kept {kept} reads and removed {removed} reads matching "
        f"{left_anchor}[{bp_between_anchors}bp]{right_anchor}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
