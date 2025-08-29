# Project Files Overview

This repository contains datasets, genome files, motif scans, and configuration files for analyzing **binding sites within the *Cni-xol-1* promoter** after knock-in (KI) into the *C. briggsae* genome. ATAC-seq analysis revealed a **specific peak in the *Cni-xol-1* promoter region**, and this repository provides motif and binding site analysis of that peak.

👉 You can **visualize all data files interactively in JBrowse2** using the following link:
[link](https://wormbase.org/tools/genome/jbrowse2/?config=https://raw.githubusercontent.com/lybCNU/xol1RI/refs/heads/main/ATACseq/config.json)

---

## File Contents

### Motif and Binding Site Files

* **tra-1 motif**
  Motif predicted using **MEME-ChIP** based on 184 *C. elegans* TRA-1 ChIP-seq sites (Berkseth et al., 2013).

* **TRA-1\_site\_PMC2778738.bed**
  TRA-1 binding sites derived from the study referenced in PMC2778738.

* **TRA-1\_site\_fimo.gff**
  Results from scanning the *cb4CnXol* genome with the **MEME-ChIP-derived TRA-1 motif** using FIMO.

* **binding\_site\_by\_jaspar.gff3**
  Binding-site predictions from JASPAR scans of the *Cni-xol-1* promoter-specific ATAC-seq peak:

  * **SEX-1 sites**: relative profile score threshold **80%**
  * **T-box sites**: relative profile score threshold **85%**

### Genome Files

* **cb4CnXol.fa.gz** / **cb4CnXol.fa.gz.gzi**
  Reference genome sequence for *C. briggsae* with *Cni-xol-1* knock-in, and its index file.

* **cb4CnXol.gff.gz** / **cb4CnXol.gff.gz.tbi**
  Genome annotation file and corresponding index.

### ATAC-seq Data

* **cb4CnXol\_mixedSex\_embryo\_ATACseq.bw**
  BigWig track of ATAC-seq signal from mixed-sex embryo samples.

* **cb4CnXol\_narrowPeak\_macs3.bed**
  NarrowPeak file of ATAC-seq peaks called by MACS3.

### Configuration and Interval Files

* **config.json**
  Configuration file used to generate JBrowse2 sessions and links.

* **thuB129.bed**
  Genomic intervals of **thuB129** aligned to the *cb4CnXol* genome.

---

## Notes

This collection of files supports visualization and analysis of **binding sites in the *Cni-xol-1* promoter knock-in locus within the *C. briggsae* genome**—integrating motif scans (TRA-1 / SEX-1 / T-box), ATAC-seq signal/peaks, and published binding sites. All datasets can be explored using the **interactive JBrowse2 link** provided above.
