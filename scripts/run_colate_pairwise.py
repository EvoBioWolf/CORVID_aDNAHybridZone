#!/usr/bin/env python3

import csv
import itertools
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime
import sys

TYP = sys.argv[2] if len(sys.argv) > 2 else "wgs"
CHR = sys.argv[1] if len(sys.argv) > 1 else "all_chr"
N_WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 1
DRY_RUN = False

# ============================================================
# SETTINGS
# ============================================================

# Folder containing the .colate.in files and chromosome file.
NUM_BOOTSTRAPS = 20
YEARS_PER_GEN = 5.79
BINS = "3,7,0.2"

# ============================================================
# PATHS
# ============================================================

dat = Path("PATH")

# Main Colate directory
colate_dir = dat / "02_results" / "colate"

# Colate executable
colate_bin = (
    dat
    / "Colate"
    / "binaries"
    / "v0.1.5_x86_64_dynamic"
    / "bin"
    / "Colate"
)

# Input .colate.in files:
# e.g. 02_results/colate/all_chr/out_E01_1x_all_chr.colate.in
input_dir = colate_dir / TYP / CHR

# Metadata table:
# 02_results/colate/heatmap/colate_samples.txt
sample_file = colate_dir / "heatmap" / TYP /"colate_samples.txt"

# Pairwise .coal outputs go here:
# 02_results/colate/heatmap/
output_dir = colate_dir / "heatmap" / TYP

# Log file
log_file = output_dir / f"colate_pairwise_{CHR}.log"

# Colate chromosome/region file:
# e.g. 02_results/colate/all_chr.txt
chr_file = colate_dir / f"{CHR}.txt"


# ============================================================
# CHECK PATHS AND CREATE OUTPUT DIRECTORY
# ============================================================

output_dir.mkdir(parents=True, exist_ok=True)

if not colate_bin.exists():
    raise FileNotFoundError(
        f"Colate binary was not found:\n{colate_bin}"
    )

if not sample_file.exists():
    raise FileNotFoundError(
        f"Sample metadata file was not found:\n{sample_file}"
    )

if not input_dir.exists():
    raise FileNotFoundError(
        f"Input directory was not found:\n{input_dir}"
    )

if not chr_file.exists():
    raise FileNotFoundError(
        f"Chromosome/region file was not found:\n{chr_file}"
    )


# ============================================================
# READ SAMPLE METADATA
# ============================================================

samples = []

with open(sample_file, newline="") as handle:

    reader = csv.DictReader(handle, delimiter="\t")

    required_columns = {"sample", "age_years", "population"}

    if reader.fieldnames is None:
        raise ValueError(
            f"Metadata file has no header:\n{sample_file}"
        )

    if not required_columns.issubset(reader.fieldnames):
        raise ValueError(
            f"Expected columns: {required_columns}\n"
            f"Found columns: {reader.fieldnames}"
        )

    for row in reader:

        sample_id = row["sample"].strip()
        age_years = row["age_years"].strip()
        population = row["population"].strip()

        if not sample_id:
            continue

        try:
            float(age_years)
        except ValueError:
            raise ValueError(
                f"Invalid age_years value for {sample_id}: {age_years}"
            )

        samples.append(
            {
                "sample": sample_id,
                "age_years": age_years,
                "population": population,
            }
        )

if len(samples) < 2:
    raise ValueError(
        f"Need at least two samples; found {len(samples)}."
    )

print(f"Read {len(samples)} samples from:")
print(sample_file)

# Unique unordered sample pairs:
# A-vs-B is included once; B-vs-A is not separately included.
pairs = list(itertools.combinations(samples, 2))

print(f"Number of unique pairs: {len(pairs)}")
print(f"Input group: {CHR}")
print(f"Workers: {N_WORKERS}")
print(f"Dry run: {DRY_RUN}")


# ============================================================
# CHECK THAT ALL .COLATE.IN FILES EXIST
# ============================================================

missing_tmp = []

for sample_info in samples:

    sample_id = sample_info["sample"]

    tmp_file = input_dir / f"out_{sample_id}_{CHR}.colate.in"

    if not tmp_file.exists() or tmp_file.stat().st_size == 0:
        missing_tmp.append(tmp_file)

if missing_tmp:

    print("\nERROR: Missing or empty .colate.in files.")

    for tmp_file in missing_tmp[:20]:
        print(tmp_file)

    if len(missing_tmp) > 20:
        print(f"... plus {len(missing_tmp) - 20} more missing files.")

    raise SystemExit(1)

print("All required .colate.in files were found.")


# ============================================================
# RUN ONE PAIRWISE COLATE ANALYSIS
# ============================================================

def run_pair(pair):

    target, reference = pair

    target_id = target["sample"]
    target_age = target["age_years"]

    reference_id = reference["sample"]
    reference_age = reference["age_years"]

    # Input files remain in e.g. colate/all_chr/
    target_tmp = input_dir / f"out_{target_id}_{CHR}.colate.in"
    reference_tmp = input_dir / f"out_{reference_id}_{CHR}.colate.in"

    # Outputs go into colate/heatmap/
    output_prefix = output_dir / (
        f"pairwise_{target_id}_vs_{reference_id}_{CHR}"
    )

    output_coal = Path(str(output_prefix) + ".coal")

    # Restart-safe: skip existing non-empty output files.
    if output_coal.exists() and output_coal.stat().st_size > 0:
        return {
            "status": "SKIPPED",
            "target": target_id,
            "reference": reference_id,
            "message": str(output_coal),
        }

    command = [
        str(colate_bin),
        "--mode", "mut",
        "--mut", "./data/mut-ages/relate",
        "--target_tmp", str(target_tmp),
        "--reference_tmp", str(reference_tmp),
        "--bins", BINS,
        "--chr", f"{CHR}.txt",
        "--num_bootstraps", str(NUM_BOOTSTRAPS),
        "--target_age", str(target_age),
        "--reference_age", str(reference_age),
        "--years_per_gen", str(YEARS_PER_GEN),
        "-o", str(output_prefix),
    ]

    command_text = " ".join(command)

    if DRY_RUN:
        return {
            "status": "DRY_RUN",
            "target": target_id,
            "reference": reference_id,
            "message": command_text,
        }

    try:

        result = subprocess.run(
            command,
            cwd=colate_dir,
            capture_output=True,
            text=True,
            check=True,
        )

        # Confirm an output was actually created.
        if not output_coal.exists() or output_coal.stat().st_size == 0:
            return {
                "status": "FAILED",
                "target": target_id,
                "reference": reference_id,
                "message": (
                    "Colate returned success but .coal output is missing/empty.\n"
                    f"Expected: {output_coal}\n"
                    f"STDOUT:\n{result.stdout}\n"
                    f"STDERR:\n{result.stderr}"
                ),
            }

        return {
            "status": "DONE",
            "target": target_id,
            "reference": reference_id,
            "message": str(output_coal),
        }

    except subprocess.CalledProcessError as error:

        return {
            "status": "FAILED",
            "target": target_id,
            "reference": reference_id,
            "message": (
                f"Command:\n{command_text}\n\n"
                f"STDOUT:\n{error.stdout}\n\n"
                f"STDERR:\n{error.stderr}"
            ),
        }


# ============================================================
# RUN ALL PAIRS
# ============================================================

start_time = datetime.now()

n_done = 0
n_skipped = 0
n_failed = 0

with open(log_file, "a") as log:

    log.write("\n")
    log.write("=" * 80 + "\n")
    log.write(f"Started: {start_time}\n")
    log.write(f"CHR/input group: {CHR}\n")
    log.write(f"Samples: {len(samples)}\n")
    log.write(f"Pairs: {len(pairs)}\n")
    log.write(f"Workers: {N_WORKERS}\n")
    log.write(f"Dry run: {DRY_RUN}\n")
    log.write("=" * 80 + "\n")

    with ThreadPoolExecutor(max_workers=N_WORKERS) as executor:

        futures = {
            executor.submit(run_pair, pair): pair
            for pair in pairs
        }

        for future in as_completed(futures):

            result = future.result()

            status = result["status"]
            target = result["target"]
            reference = result["reference"]
            message = result["message"]

            print(f"{status}: {target} vs {reference}")

            # Keep log one pair at a time; replace line breaks for easy reading.
            log_message = message.replace("\n", " | ")

            log.write(
                f"{status}\t{target}\t{reference}\t{log_message}\n"
            )

            log.flush()

            if status == "DONE":
                n_done += 1

            elif status == "SKIPPED":
                n_skipped += 1

            elif status == "FAILED":
                n_failed += 1


# ============================================================
# FINAL SUMMARY
# ============================================================

end_time = datetime.now()

print("\n" + "=" * 80)
print("Pairwise Colate analysis finished.")
print(f"Start time: {start_time}")
print(f"End time: {end_time}")
print(f"Completed: {n_done}")
print(f"Skipped existing outputs: {n_skipped}")
print(f"Failed: {n_failed}")
print(f"Log file: {log_file}")
print("=" * 80)

if n_failed > 0:
    raise SystemExit(1)
