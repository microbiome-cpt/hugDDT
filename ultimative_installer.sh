#!/usr/bin/bash
set -euo pipefail

# --- init conda ---
if ! command -v conda &>/dev/null; then
    echo "conda not found in PATH" >&2
    exit 1
fi

# ensure conda shell functions
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "LETS INSTALL CONDA ENVS!"
# --- create environments ---
for yml in envs/*.yml; do
    mamba env create -f "$yml"
done

echo "CHECK IF THE ENVS ARE IN THE LIST:"
conda env list
echo ""

# --- install third-party tools ---
echo "LETS INSTALL DORADO!"
if ! command -v dorado &>/dev/null; then
    wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.2.0-linux-x64.tar.gz
    tar xzf dorado-1.2.0-linux-x64.tar.gz
    echo "export PATH=$(pwd)/dorado-1.2.0-linux-x64/bin:$PATH" >> "$CONDA_PREFIX/etc/conda/activate.d/dorado.sh"
	export PATH=$(pwd)/dorado-1.2.0-linux-x64/bin:$PATH
fi
echo ""
echo "CHECK IF DORADO VERSION APPEARS:"
dorado -v
echo ""

echo "LETS INSTALL NEXTPOLISH!"
# --- NextPolish ---
if ! command -v nextpolish &>/dev/null; then
    wget https://github.com/Nextomics/NextPolish/releases/download/v1.4.1/NextPolish.tgz
    tar -vxzf NextPolish.tgz
    echo "export PATH=\$PATH:$(pwd)/NextPolish" >> "$CONDA_PREFIX/etc/conda/activate.d/nextpolish.sh"
    cd NextPolish && make && cd ../
	export PATH=$(pwd)/NextPolish:$PATH
fi

echo "DON'T BOTHER IF YOU SEE SOME [samtools] Make Error.."
echo ""
echo "CHECK IF NEXTPOLISH VERSION APPAERS:"
nextPolish -v
