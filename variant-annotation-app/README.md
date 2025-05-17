# Variant Annotation App

A Python-based variant annotation application integrating multiple popular annotators: **bcftools**, **snpEff**, **Ensembl VEP**, and **ANNOVAR**.

---

## 🚀 Features

- Annotate VCF files using **bcftools annotate**.
- Annotate variants with **snpEff** (Java-based).
- Annotate variants with **Ensembl Variant Effect Predictor (VEP)** (Perl-based).
- Annotate variants with **ANNOVAR** (command-line Perl tool).
- Parse annotated VCF outputs into Python dictionaries for downstream processing.

---

## 📦 Prerequisites

- Python 3.x
- **bcftools** installed and accessible in your system PATH
- Java Runtime Environment (JRE) for running snpEff
- Perl (v5.10 or higher) with required Ensembl Perl API modules for VEP
- **ANNOVAR** downloaded and configured
- Git (to clone Ensembl API modules if needed)

---

## ⚙️ Setup Instructions

### 1. Install bcftools

```bash
# On macOS using Homebrew:
brew install bcftools
```

---

### 2. Download snpEff

- Download `snpEff.jar` and place it in:

```
/path/to/variant-annotation-app/bin/snpEff/snpEff.jar
```

- Ensure Java is installed:

```bash
java -version
```

---

### 3. Download and Setup Ensembl VEP

- Download Ensembl VEP from the [official site](https://www.ensembl.org/info/docs/tools/vep/index.html) and extract it to:

```
/path/to/variant-annotation-app/bin/ensembl-vep-release-114/
```

- Clone the required Ensembl Perl API repositories:

```bash
cd /path/to/variant-annotation-app/bin/
git clone https://github.com/Ensembl/ensembl.git
git clone https://github.com/Ensembl/ensembl-variation.git
git clone https://github.com/Ensembl/ensembl-funcgen.git
git clone https://github.com/Ensembl/ensembl-external.git
```

- Set the `PERL5LIB` environment variable:

```bash
export PERL5LIB=/path/to/variant-annotation-app/bin/ensembl/modules:\
/path/to/variant-annotation-app/bin/ensembl-variation/modules:\
/path/to/variant-annotation-app/bin/ensembl-funcgen/modules:\
/path/to/variant-annotation-app/bin/ensembl-external/modules:\
/path/to/variant-annotation-app/bin/ensembl-vep-release-114/modules:$PERL5LIB
```

- Add the `export` line to your `.bashrc`, `.bash_profile`, or `.zshrc` for persistence.

---

### 4. Download and Setup ANNOVAR

- Download ANNOVAR from [http://annovar.openbioinformatics.org](http://annovar.openbioinformatics.org).

- Extract it to:

```
/path/to/variant-annotation-app/bin/annovar/
```

- Confirm Perl is installed:

```bash
perl -v
```

---

## 🐳 Docker Setup

You can also build and run the app using Docker:

```bash
docker build -t variant-annotation-app .
docker run -it -p 8000:8000 variant-annotation-app
```

Then access your app at: [http://localhost:8000](http://localhost:8000)

---

## 🛠️ Troubleshooting

- **Missing Perl modules for VEP**: Clone Ensembl API repos and correctly set `PERL5LIB`.
- **Java not found for snpEff**: Ensure Java JRE is installed and added to PATH.
- **bcftools not found**: Install via Homebrew or package manager and confirm it's in PATH.
- **ANNOVAR issues**: Confirm Perl is installed, ANNOVAR scripts are executable, and database files are downloaded properly.

---
