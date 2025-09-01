# Cni-xol-1 Gene Model Revision

This repository contains files related to the transcript structure reanalysis of *Caenorhabditis nigoni* genes **Cni-xol-1** and **Cni-xol-1.1**.

## Contents

* **Public\_RNAseq\_on\_Cni-xol-1.bam / .bam.bai**: RNA-seq alignments.
* **Cni-xol-1\_annotation\_update.gff3**: Updated gene model from StringTie + TransDecoder.
* **Cni-xol-1vsCni-xol-1.1.paf**: Pairwise alignment of *Cni-xol-1* vs *Cni-xol-1.1*.
* **cb4\_vs\_JU1421\_ASM2792064v1.paf**: Genome alignment for visualization.

## Visualization

Interactive genome browser: [View in JBrowse2](https://wormbase.org/tools/genome/jbrowse2/?config=https%3A%2F%2Fraw.githubusercontent.com%2FlybCNU%2Fxol1RI%2Frefs%2Fheads%2Fmain%2FCni-xol-1_GeneModel_Revision%2Fconfig.json)

## Methods

* RNA-seq aligned with **HISAT2 v2.2.1**.
* BAM processed with **SAMtools v1.21**.
* Transcripts assembled with **StringTie v2.2.3**.
* CDS predicted with **TransDecoder v5.7.1**.
* Alignments generated with **minimap2 (asm20)**.

---

For details, see associated analyses and manuscript.
