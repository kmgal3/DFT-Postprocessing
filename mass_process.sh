#!/bin/bash
# Usage: ./mass_process.sh max_jobs parent_folder calc_folder [args...]
# Runs ./vaspout_h5.py inside every matching subdir, keeping at most max_jobs in parallel.

MAX_JOBS=$1
parent_folder=$2
calc_folder=$3
shift 3
script_args=("$@")

if [[ -z "$MAX_JOBS" || -z "$parent_folder" || -z "$calc_folder" ]]; then
    echo "Usage: $0 max_jobs parent_folder calc_folder [args...]"
    exit 1
fi

# Fixed script to execute inside each subdir
START_DIR="$(pwd)"
SCRIPT="$START_DIR/vaspout_h5.py"

# Function to throttle jobs to MAX_JOBS
function throttle() {
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        sleep 1
    done
}

# Find all subdirectories named $calc_folder and run in parallel
find "$parent_folder" -type d -name "$calc_folder" | while read -r dir; do
    throttle
    (
        echo ">>> Entering $dir"
        cd "$dir" || exit 1
        "python3" "$SCRIPT" "${script_args[@]}"
    ) &
done

wait
echo "All jobs finished."
