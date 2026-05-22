import os

# --- Configuration ---
run = "399465"
output_py = f"filelists/files_{run}_cfi.py"  # output file name

base_dir = "/eos/cms/store/t0streamer/Data/"
stream_prefix = "PhysicsHITrackerNZS"      # NZS0, NZS1, ..., NZS9
file_extension = ".dat"

# Build the subdirectories descending from run number
subdir = os.path.join("000", run[:3], run[3:])

# --- Gather all input dirs (NZS0 ... NZS9) ---
input_dirs = [
    os.path.join(base_dir, f"{stream_prefix}{i}", subdir)
    for i in range(10)
]

all_files = []

# --- Loop over all NZS directories ---
for d in input_dirs:
    abs_dir = os.path.abspath(d)
    if not os.path.exists(abs_dir):
        print(f"Directory not found: {abs_dir}")
        continue

    files = sorted([
        os.path.join(abs_dir, f)
        for f in os.listdir(abs_dir)
        if os.path.isfile(os.path.join(abs_dir, f)) and f.endswith(file_extension)
    ])

    if not files:
        print(f"No {file_extension} files found in {abs_dir}")
    else:
        print(f"Found {len(files)} files in {abs_dir}")
        all_files.extend(files)

if not all_files:
    raise RuntimeError("No input files found in any NZS directory.")

# --- Write output python file ---
with open(output_py, "w") as f:
    f.write("import FWCore.ParameterSet.Config as cms\n\n")
    f.write("readFiles = cms.untracked.vstring(*(\n")

    for full_path in all_files:
        f.write(f"    'file:{full_path}',\n")

    f.write("))\n")

print(f"Wrote {len(all_files)} files to {output_py}")
