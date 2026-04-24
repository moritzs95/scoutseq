"""Load downstream DRO target parameters from pipeline-exported environment variables."""

import os
from typing import Optional


def _require_env(name: str) -> str:
    value = os.environ.get(name, "").strip()
    if not value:
        raise SystemExit(f"Required DRO config variable '{name}' is missing.")
    return value


def _int_env(name: str, default: Optional[int] = None) -> int:
    value = os.environ.get(name, "").strip()
    if not value:
        if default is None:
            raise SystemExit(f"Required DRO config variable '{name}' is missing.")
        return default
    try:
        return int(value)
    except ValueError as exc:
        raise SystemExit(f"DRO config variable '{name}' must be an integer, got '{value}'.") from exc


def load_dro_config(*, require_gs_primer: bool = False, require_lengths: bool = False, require_insertion_type: bool = False) -> dict[str, object]:
    hdr_anchor = os.environ.get("SCOUT_DRO_HDR_ANCHOR", "").strip()
    insertion_type = os.environ.get("SCOUT_DRO_INSERTION_TYPE", "substitution").strip() or "substitution"
    config: dict[str, object] = {
        "label": os.environ.get("SCOUT_DRO_LABEL", "").strip(),
        "left_anchor": _require_env("SCOUT_DRO_LEFT_ANCHOR"),
        "right_anchor": _require_env("SCOUT_DRO_RIGHT_ANCHOR"),
        "hdr_anchor": hdr_anchor,
        "wt_amplicon": _require_env("SCOUT_DRO_WT_AMPLICON"),
        "wt_amplicon_after_cutsite": os.environ.get("SCOUT_DRO_WT_AMPLICON_AFTER_CUTSITE", "").strip()
        or _require_env("SCOUT_DRO_WT_AMPLICON"),
        "hdrbc_len": _int_env("SCOUT_DRO_HDRBC_LENGTH"),
        "thresholds": {
            "I90plus": _int_env("SCOUT_DRO_THRESHOLD_I90PLUS", 30),
            "D90plus": _int_env("SCOUT_DRO_THRESHOLD_D90PLUS", 30),
            "S1plus": _int_env("SCOUT_DRO_THRESHOLD_S1PLUS", 50),
        },
    }

    if require_gs_primer:
        gs_primer = _require_env("SCOUT_DRO_GS_PRIMER")
        config["gs_primer"] = gs_primer
        config["gs_primer_length"] = len(gs_primer)

    if require_lengths:
        config["cdna_length"] = _int_env("SCOUT_DRO_CDNA_LENGTH")
        config["reference_length"] = _int_env("SCOUT_DRO_REFERENCE_LENGTH")

    if require_insertion_type:
        if insertion_type not in {"substitution", "cINS", "rINS"}:
            raise SystemExit(
                "DRO config variable 'SCOUT_DRO_INSERTION_TYPE' must be one of: substitution, cINS, rINS."
            )
        if not hdr_anchor and insertion_type != "rINS":
            raise SystemExit("Required DRO config variable 'SCOUT_DRO_HDR_ANCHOR' is missing.")
        config["insertion_type"] = insertion_type
        if insertion_type == "cINS":
            config["cins_hdrbc"] = _require_env("SCOUT_DRO_CINS_HDRBC")
    elif not hdr_anchor and insertion_type != "rINS":
        raise SystemExit("Required DRO config variable 'SCOUT_DRO_HDR_ANCHOR' is missing.")

    return config
