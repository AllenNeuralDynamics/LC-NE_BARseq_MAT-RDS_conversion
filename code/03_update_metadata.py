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

import shutil
from datetime import datetime, timezone
from pathlib import Path

from aind_data_schema.components.identifiers import Code
from aind_data_schema.core.data_description import DataDescription
from aind_data_schema.core.processing import DataProcess, ProcessStage, Processing
from aind_data_schema_models.process_names import ProcessName

PROCESS_NAME_LABEL = "processed_MAT2RDS"
RESULTS_DIR = Path("/results")
DATA_DIR = Path("/data")
EXPERIMENTERS = ["Polina Kosillo"]

# Code Ocean capsule URL recorded in processing.json. Hardcoded to the released
# capsule so the published metadata always references the canonical capsule,
# regardless of which working/test capsule is running this script.
CAPSULE_URL = "https://codeocean.allenneuraldynamics.org/capsule/3953531/tree"


def find_output_dirs() -> list[Path]:
    """Return the per-subject output directories created by the conversion scripts."""
    return sorted(
        p for p in RESULTS_DIR.iterdir()
        if p.is_dir() and "_processed_MAT2RDS_" in p.name
    )


def parse_input_asset_name(out_dir_name: str) -> str:
    """Strip the '_processed_MAT2RDS_<timestamp>' suffix to recover the input asset name."""
    return out_dir_name.split("_processed_MAT2RDS_")[0]


def copy_peer_metadata(input_dir: Path, out_dir: Path) -> None:
    """Copy acquisition/procedures/subject JSON files from the input asset, unchanged.

    Raises FileNotFoundError if any expected file is missing.
    """
    for fname in ("acquisition.json", "procedures.json", "subject.json"):
        src = input_dir / fname
        if not src.exists():
            raise FileNotFoundError(f"Required peer metadata file not found: {src}")
        shutil.copy(src, out_dir / fname)
        print(f"  copied {fname}")


def write_data_description(input_dir: Path, out_dir: Path) -> None:
    """Build and write a derived data_description.json from the input asset's raw one.

    Uses ``DataDescription.from_data_description`` so the schema sets
    ``data_level=DERIVED``, populates ``source_data`` with the input asset
    name, and generates an AIND-compliant derived asset name.
    """
    raw_dd_path = input_dir / "data_description.json"
    raw_dd = DataDescription.model_validate_json(raw_dd_path.read_text())
    derived_dd = DataDescription.from_data_description(
        raw_dd, process_name=PROCESS_NAME_LABEL
    )
    derived_dd.write_standard_file(output_directory=out_dir)
    print(f"  wrote derived data_description.json (source: {raw_dd.name})")


def write_processing(out_dir: Path) -> None:
    """Build and write a processing.json describing this conversion run."""
    code = Code(
        url=CAPSULE_URL,
        name="LC-NE_BARseq_MAT-RDS_conversion",
        run_script=Path("code/run"),
        language="R",
    )
    process = DataProcess(
        process_type=ProcessName.FILE_FORMAT_CONVERSION,
        name="MAT to RDS conversion",
        stage=ProcessStage.PROCESSING,
        code=code,
        experimenters=EXPERIMENTERS,
        start_date_time=datetime.now(timezone.utc),
        notes="Converts BARseq .mat (HDF5 v7.3) to SingleCellExperiment .rds files.",
    )
    processing = Processing(data_processes=[process])
    processing.write_standard_file(output_directory=out_dir)
    print("  wrote processing.json")


def update_metadata_for_subject(out_dir: Path) -> None:
    """Generate AIND-compliant metadata files for one per-subject output folder."""
    input_asset = parse_input_asset_name(out_dir.name)
    input_dir = DATA_DIR / input_asset

    print(f"Updating metadata for {out_dir.name}")
    print(f"  input asset: {input_asset}")

    copy_peer_metadata(input_dir, out_dir)
    write_data_description(input_dir, out_dir)
    write_processing(out_dir)


def main() -> None:
    """Discover output folders under /results/ and generate metadata for each."""
    print(f"Capsule URL: {CAPSULE_URL}\n")

    out_dirs = find_output_dirs()
    if not out_dirs:
        raise RuntimeError(f"No output directories found under {RESULTS_DIR}")
    print(f"Found {len(out_dirs)} output director(ies)")
    for out_dir in out_dirs:
        print()
        update_metadata_for_subject(out_dir)
    print()
    print("Done.")


if __name__ == "__main__":
    main()
