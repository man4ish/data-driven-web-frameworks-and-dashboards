#!/bin/bash

# Directory where SnpEff data should go
SNPEFF_DIR="/opt/snpEff/data"

mkdir -p "$SNPEFF_DIR"

echo "Downloading SnpEff hg38 database..."

wget -qO- "https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_GRCh38.99.zip" \
    | unzip -o -d "$SNPEFF_DIR" -

echo "SnpEff hg38 database downloaded to $SNPEFF_DIR"

