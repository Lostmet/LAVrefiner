# LAVrefiner (SVrefiner)

**LAVrefiner** is a Python-based tool for processing structural variants (SVs) and indels (Length-altering variants, LAVs) and generating refined LAVs (rLAVs) in VCF v4.2 format.
It integrates sequence alignment, window scanning, and merging workflows to resolve overlapping LAVs, with a focus on `DEL` (deletions) and `INS` (insertions). It also supports Tandem Repeat (TR) aware alignments, Quality Control (QC), and TR motif dosage calculations.

Complex LAVs (i.e., where ref and alt differ at the first base, or both ref and alt are longer than 1 bp) are not processed in the current version to improve efficiency.

> ref: reference allele, alt: alternative allele

---

## Quick Start

### Installation (Conda)

```bash
# Create and activate environment
conda create -n LAVrefiner python=3.12
conda activate LAVrefiner

# Clone and install LAVrefiner
git clone https://github.com/Lostmet/LAVrefiner.git
cd LAVrefiner
pip install .

```

You can replace the above steps by downloading the release package from the Releases page, which already contains the pre-built executable.

```bash
# Install MAFFT
conda install conda-forge::mafft

```

Alternatively, you can download MAFFT directly from the [official MAFFT website](https://mafft.cbrc.jp/alignment/software/) and install it manually. In this case, please make sure that the `mafft` executable is accessible through your environment `PATH`.

---

## Input

**Required input files:**

* **VCF** file with index (`.vcf.gz` and `.vcf.gz.tbi` or `.vcf.gz.csi`)
* **FASTA** reference genome (`.fa` or `.fasta`) including the chromosomes covered by the LAVs.
* **(Optional) TR VCF Database**: A Tandem Repeat reference database for TR-aware processing and motif extraction (`--tr-vcf`). This can be downloaded from [TRCompDB (https://zenodo.org/records/13263615)](https://zenodo.org/records/13263615). Please note that the input must be the VCF format file from TRCompDB, and you need to select the file corresponding to your reference sequence. It is highly recommended to use the `oriMotifs` file to ensure a broader selection of Motifs.

---

## Main Output

| File/Directory | Description |
| --- | --- |
| `alignment_results/` | Alignment results for each overlapping group |
| `matrix_results/` | Output matrices for each overlapping group |
| `rLAV.vcf` | Refined LAVs (final VCF) |
| `nLAV.vcf` | Non-overlapping LAVs |
| `oLAV.vcf` | Overlapping LAVs |
| `*.log` | Execution logs |

---

### Notes

* **Chromosome identifiers in VCF and FASTA must be numeric.**
* Examples: `1`, `2`, `3` (autosomes), `23` (X), `24` (Y)
* Identifiers such as `chr1`, `chrX` are not supported.


* **LAVs must be normalized and left-aligned in standard VCF format.**
* Only **biallelic LAVs** are supported; multiallelic records are not accepted.
* If your VCF is not normalized, please preprocess it using external tools (e.g., `bcftools norm`) before running LAVrefiner.


* **Since the quality of input LAVs (especially the genotyping rate) directly affects the genotyping rate of the final rLAVs, it is strongly recommended to perform quality control (QC) before running LAVrefiner.**

---

## VCF Example

> The VCF file must be compressed (`.vcf.gz`) and indexed (`.tbi` or `.csi`).

```vcf
##fileformat=VCFv4.2
##source=YourTool
#CHROM  POS  ID    REF    ALT     QUAL    FILTER  INFO    FORMAT Sample1  Sample2  Sample3  Sample4
1       1    lav1  ACTA   A       50      PASS    .       GT      1/1      1/0      0/0      ./.
1       5    lav2  G      GAAC    99      PASS    .       GT      0/0      1/0      0/0      0/0
1       6    lav3  GCTAG  <INV>   98      PASS    .       GT      ./.      0/0      1/1      1/1

```

---

## Commands

| Command | Function |
| --- | --- |
| `process-vcf` | Process VCF to identify overlapping LAVs; outputs `nLAV.vcf`, `oLAV.vcf` and `variants_pre_aligned.fasta`. |
| `align` | Perform sequence alignment for overlapping LAVs (Hybrid MAFFT + TR-aware mode). |
| `make-rlav` | Define refined LAVs from aligned sequences and output `rLAV.vcf`. |
| `run-all` | Execute the complete pipeline. |
| `QC` | Quality Control and filtering for the final `rLAV.vcf`. |
| `dosage` | Calculate motif dosage and copy number variations for TR loci suitable for GWAS/eQTL. |

## Common Options

| Option | Description |
| --- | --- |
| `--vcf` | **Required.** Input VCF file (`.vcf.gz`; index `.csi` or `.tbi` required). |
| `--ref` | **Required.** Input reference FASTA file (`.fasta` or `.fa`). |
| `--out` | **Required.** Output directory. |
| `--threads` | *Optional.* Number of threads to use (default: `10`). |
| `--write-matrix` | *Optional.* `YES`/`NO` (default: `NO`). Whether to output **X** and **T** matrices. |
| `--tr-vcf` | *Optional.* Path to TRCompDB VCF for TR-aware processing and motif extraction. |

### Example Workflow

#### Step 1: Process VCF

```bash
lavrefiner process-vcf --vcf test.vcf.gz --ref test.fasta --out test_out --threads 10

```

This separates non-overlapping (`nLAV.vcf`) and overlapping (`oLAV.vcf`) variants and generates pre-aligned fasta sequences.

#### Step 2: Align

```bash
lavrefiner align --vcf test.vcf.gz --ref test.fasta --out test_out --threads 10 --tr-vcf TRCompDB.vcf.gz

```

Performs hybrid alignment and saves TR metadata (if `--tr-vcf` is provided).

#### Step 3: Define rLAV

```bash
lavrefiner make-rlav --vcf test.vcf.gz --ref test.fasta --out test_out --threads 10

```

Generates the final `rLAV.vcf` and associated metadata.

#### Complete Pipeline (Run All)

```bash
lavrefiner run-all --vcf test.vcf.gz --ref test.fasta --out test_out --threads 10 --tr-vcf TRCompDB.vcf.gz

```

#### Step 4: Quality Control (QC)

Filter the resulting rLAVs based on tandem repeat divergence or repeat length.

```bash
lavrefiner QC --vcf test_out/rLAV.vcf --out test_out --prefix filtered_rLAV --tr-div 0.1 --min-rl 2

```

**QC Tag Definitions:**

* `--vcf`: The path to the input rLAV VCF file (supports `.vcf.gz`).
* `--out`: The directory where the filtered output files will be saved.
* `--prefix`: The prefix name for the newly generated filtered VCF file.
* `--tr-div`: Tandem Repeat Divergence threshold. Excludes variants with a `TR_DIV >=` this specified value (e.g., `0.1`).
* `--exclude`: Excludes a specific type of TR variant. Options are `STR`, `VNTR`, or `TR` (selecting `TR` excludes both STRs and VNTRs).
* `--extract`: Extracts and keeps ONLY a specific type of TR variant. Options are `STR`, `VNTR`, or `TR`.
* `--min-rl`: The minimum Repeat Length (RL) required to keep a variant.
* `--max-rl`: The maximum Repeat Length (RL) allowed to keep a variant.
* `--tr-vcf`: *(Optional)* Path to the TRCompDB VCF. Used for online annotation and filtering of variants that lack `TR_TYPE` or `RL` information in the original file.

#### Step 5: Dosage Calculation

Calculate TR motif dosages and copy numbers directly for downstream applications like GWAS or eQTLs. The inputed VCF should incorporate all indels and SVs (LAVs) within TR regions.

```bash
lavrefiner dosage --vcf test_out/LAV.vcf.gz --tr-vcf TRCompDB.vcf.gz --out test_out --prefix TR_dosage --threads 10

```

---

## Requirements

* **Python ≥ 3.8**
* **MAFFT ≥ v7.526**
* Python libraries:
* `pandas`
* `numpy`
* `biopython`
* `pysam`
* `click`
* `tqdm`
* `Levenshtein`



---

## Conceptual Overview

### Main
The pipeline resolves overlapping LAVs into more precise rLAVs through alignment and window-based clustering.

<p align="center">
<img src="https://github.com/user-attachments/assets/8e8de218-7207-40a9-8ddc-fb5f23fb85f5" width="500">
</p>

### Module overview

<p align="center">
<img src="https://github.com/user-attachments/assets/d7def476-6db8-446f-a326-1436d7db0e40" width="800">
</p>


---

## Output Structure

### Main Directory

| File/Directory | Description |
| --- | --- |
| `rLAV.vcf` | Final refined LAVs |
| `rLAV_meta.csv` | rLAV metadata |
| `nLAV.vcf` | Non-overlapping LAVs |
| `oLAV.vcf` | Overlapping LAVs |
| `*.log` | Runtime logs |

Example log:

```
Command:  lavrefiner run-all
  --ref   test.fasta 
  --vcf   test.vcf.gz 
  --out   test_out 
Start time:                 2025-07-07 13:24:57
End time:                   2025-07-07 13:24:59
Total runtime:              0:00:01

Total variants:               100
INV count:                    1
nLAV count:                   72
Excluded SNP count:           0
Overlapping LAVs:             28
Overlapping LAVs percentage:  28.00%
Total variant groups:         9
Final rLAV count:             45
Mode:                         TR-Aware

```

### `alignment_results/`

| File | Description |
| --- | --- |
| `Group_input_origin.fasta` | Original group sequences |
| `Group_input_sliced.fasta` | Sliced insertion sequences |
| `Group_aligned_sliced.fasta` | Aligned slices |
| `Group_aligned.fasta` | Aligned sequences |
| `Group_aligned.tr` | TR layout annotation files |

### `matrix_results/`

| File | Description |
| --- | --- |
| `Group_D_matrix.csv` | D matrix (rLAV–LAV relationships) |
| `Group_X_matrix.csv` | X matrix (LAV–sample relationships) |
| `Group_T_matrix.csv` | T matrix (rLAV–sample relationships) |

### `alignment_error_logs/`

Contains MAFFT alignment error logs, if any.

---

## Support

If you encounter issues, please open an issue on GitHub or contact the authors at:

* [fenglostmet@tju.edu.cn](https://www.google.com/search?q=mailto%3Afenglostmet%40tju.edu.cn)
* [xia_xiaoxuan@outlook.com](https://www.google.com/search?q=mailto%3Axia_xiaoxuan%40outlook.com)

---

## Citation

If you use this software, please cite:

Xia X., Wu J., Gao Z. *et al.*, Modeling structural variations sequencing information to address missing heritability and enhance risk prediction (2025) *bioRxiv*. doi:[10.1101/2025.08.07.669060](https://doi.org/10.1101/2025.08.07.669060)
