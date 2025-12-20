#!/usr/bin/bash
set -euo pipefail

# --- init conda ---
if ! command -v conda &>/dev/null; then
    echo "conda not found in PATH" >&2
    exit 1
fi

# ensure conda shell functions
source "$(conda info --base)/etc/profile.d/conda.sh"

# --- create environments ---
for yml in envs/*.yml; do
    mamba env create -f "$yml"
done


# --- install third-party tools ---

# DORADO
if ! command -v dorado &>/dev/null; then
    wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.2.0-linux-x64.tar.gz
    tar xzf dorado-1.2.0-linux-x64.tar.gz
    export PATH=$(pwd)/dorado-1.2.0-linux-x64/bin:$PATH

# --- NextPolish ---
if ! command -v nextpolish &>/dev/null; then
    git clone https://github.com/Nextomics/NextPolish.git
    cd NextPolish
    make
    cd ..
    echo "export PATH=\$PATH:$(pwd)/NextPolish" >> "$CONDA_PREFIX/etc/conda/activate.d/nextpolish.sh"
fi