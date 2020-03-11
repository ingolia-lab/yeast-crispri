---
layout: default
title: Home
---

<h1 class="page-title">Yeast CRISPRi</h1>

This site is a resource for our system for genome-wide CRISPRi
screening in budding yeast (_S. cerevisiae_), described in [McGlincy
et al., 2020](https://biorxiv.org/SUBMISSION-PENDING) (bioR&chi;iv
submission pending) . This system comprises:

* A guide RNA library with ten guides per gene for most genes

* Inducible CRISPRi implemented by

  - Expression of the dCas9-Mxi effector, and the tetR regulatory
    protein, from an integrated construct (derived from [Smith _et
    al._, _Genome Biol_
    2016](https://doi.org/10.1186/s13059-016-0900-9))

  - A library of episomal plasmids for tet-inducible expression of
    guide RNAs from our library

* Random nucleotide barcodes linked with each guide RNA

  - Multiple barcodes per guide provide independent measurements of
    guide effects

  - Linear barcode amplification by in vitro transcription, reducing
    noise relative to exponential PCR

