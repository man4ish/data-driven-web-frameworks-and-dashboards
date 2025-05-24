#!/bin/bash

# Directory where ANNOVAR humandb should go
ANNOVAR_DB_DIR="/opt/annovar/humandb"

mkdir -p "$ANNOVAR_DB_DIR"

echo "Downloading ANNOVAR hg38 databases..."

cd "$ANNOVAR_DB_DIR"

# Download example ANNOVAR database files for hg38; you can add others similarly
wget -c http://www.openbioinformatics.org/annovar/download/hg38_refGene.txt.gz
wget -c http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp35a.txt.gz
wget -c http://www.openbioinformatics.org/annovar/download/hg38_gnomad211_exome.txt.gz

echo "Unzipping downloaded files..."

gunzip -f *.gz

echo "ANNOVAR hg38 databases downloaded to $ANNOVAR_DB_DIR"

