"""Generate AIND-compliant derived-asset metadata for each per-subject output folder.

For every /results/<input_asset>_processed_MAT2RDS_<timestamp>/ directory:
  - Copy peer metadata files (acquisition.json, procedures.json, subject.json)
    from the input asset under /data/<input_asset>/ unchanged.
  - Build a derived data_description.json from the input's data_description.json
    using DataDescription.from_data_description (sets data_level=DERIVED,
    populates source_data, generates the AIND-compliant derived name).
  - Write a processing.json describing this conversion run.

All schema construction goes through aind-data-schema's Pydantic models, so
any schema violation raises and aborts the run.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from aind_data_schema.components.identifiers import Code
from aind_data_schema.core.data_description import DataDescription
from aind_data_schema.core.processing import DataProcess, ProcessStage, Processing
from aind_data_schema_models.process_names import ProcessName

PROCESS_NAME_LABEL = "processed_MAT2RDS"
RESULTS_DIR = Path("/results")
DATA_DIR = Path("/data")
PEER_METADATA_FILES = ("acquisition.json", "procedures.json", "subject.json")
EXPERIMENTERS = ["Polina Kosillo"]

# Code Ocean capsule URL recorded in processing.json. Hardcoded to the released
# capsule so the published metadata always references the canonical capsule,
# regardless of which working/test capsule is running this script.
CAPSULE_URL = "https://codeocean.allenneuraldynamics.org/capsule/3953531/tree"


def get_commit_hash() -> str | None:
    """Return the capsule's git commit hash, or None if unavailable."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd="/code",
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip() or None
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def find_output_dirs() -> list[Path]:
    return sorted(
        p for p in RESULTS_DIR.iterdir()
        if p.is_dir() and "_processed_MAT2RDS_" in p.name
    )


def parse_input_asset_name(out_dir_name: str) -> str:
    """'barseq_780345_..._12-00-00_processed_MAT2RDS_<ts>' -> 'barseq_780345_..._12-00-00'"""
    return out_dir_name.split("_processed_MAT2RDS_")[0]


def update_metadata_for_subject(
    out_dir: Path, commit_hash: str | None
) -> None:
    input_asset = parse_input_asset_name(out_dir.name)
    input_dir = DATA_DIR / input_asset

    print(f"Updating metadata for {out_dir.name}")
    print(f"  input asset: {input_asset}")

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input asset directory not found: {input_dir}")

    for fname in PEER_METADATA_FILES:
        src = input_dir / fname
        if src.exists():
            shutil.copy(src, out_dir / fname)
            print(f"  copied {fname}")
        else:
            print(f"  WARN: {src} not found; skipping")

    raw_dd_path = input_dir / "data_description.json"
    if not raw_dd_path.exists():
        raise FileNotFoundError(
            f"{raw_dd_path} required to build derived data_description.json"
        )

    raw_dd = DataDescription.model_validate_json(raw_dd_path.read_text())
    derived_dd = DataDescription.from_data_description(
        raw_dd, process_name=PROCESS_NAME_LABEL
    )
    (out_dir / "data_description.json").write_text(
        derived_dd.model_dump_json(indent=2)
    )
    print(f"  wrote derived data_description.json (source: {raw_dd.name})")

    now = datetime.now(timezone.utc)
    code = Code(
        url=CAPSULE_URL,
        name="LC-NE_BARseq_MAT-RDS_conversion",
        # aind-data-schema 2.6.0's Code has no commit_hash field yet (added
        # post-2.6.0 on master); using `version` for the git commit hash.
        version=commit_hash,
        run_script=Path("code/run"),
        language="R",
    )
    process = DataProcess(
        process_type=ProcessName.FILE_FORMAT_CONVERSION,
        name="MAT to RDS conversion",
        stage=ProcessStage.PROCESSING,
        code=code,
        experimenters=EXPERIMENTERS,
        start_date_time=now,
        notes="Converts BARseq .mat (HDF5 v7.3) to SingleCellExperiment .rds files.",
    )
    processing = Processing(data_processes=[process])
    (out_dir / "processing.json").write_text(processing.model_dump_json(indent=2))
    print("  wrote processing.json")


def main() -> int:
    commit_hash = get_commit_hash()
    print(f"Capsule URL: {CAPSULE_URL}")
    print(f"Commit hash: {commit_hash or '(unavailable)'}")
    print()

    out_dirs = find_output_dirs()
    if not out_dirs:
        print(f"No output directories found under {RESULTS_DIR}", file=sys.stderr)
        return 1
    print(f"Found {len(out_dirs)} output director(ies)")
    for out_dir in out_dirs:
        print()
        update_metadata_for_subject(out_dir, commit_hash)
    print()
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
