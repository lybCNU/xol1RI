# README

## Project Overview

This project aims to investigate the genetic basis of hybrid male incompatibility in *Caenorhabditis*. Specifically, it focuses on how disrupted regulation of **xol-1**, the X-linked master switch for sex determination and dosage compensation, contributes to embryonic inviability in hybrid males. By analyzing transcriptomic data from reciprocal crosses of *C. briggsae* and *C. nigoni*, as well as comparisons with *C. latens* and *C. remanei*, the study seeks to reveal how aberrant X-chromosome repression and reduced **Cbr-xol-1** expression underlie hybrid male lethality and drive the evolution of zygotic reproductive isolation.

## Data

* **RNAseq.RDS**: Processed RNA-seq dataset containing transcriptome profiles of *C. briggsae*, *C. nigoni*, and their hybrid embryos.
* **Data\_S1.txt**: Table S1, listing *C. briggsae* and *C. nigoni* genes whose RNA-seq reads have less than 1% interspecific mapping contamination.

## Code

* **fig2\_figS2.R**: Script for generating Figure 2 and Supplementary Figure S2.
* **fig3CDE\_figS3.R**: Script for generating Figure 3CDE and Supplementary Figure S3.

## Reference Files

* **mapping\_reference/**: Directory containing reference genomes and GTF files used for sequence alignment.

## Usage

1. Load the processed RNA-seq dataset (`RNAseq.RDS`) into R.
2. Use the scripts (`fig2_figS2.R`, `fig3CDE_figS3.R`) to reproduce the figures presented in the study.
3. Consult `Data_S1.txt` for the list of genes whose RNA-seq reads have less than 1% interspecific mapping contamination.

## Citation

If you use this repository or its data in your research, please cite the corresponding study.

---

Maintainer: \[Your Name]
Contact: \[Your Email]
