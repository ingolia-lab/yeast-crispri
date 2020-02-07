---
layout: page
title: Guides v1
---

The v1 guide RNA library was designed according to empirical data
regarding effective guides for yeast CRISPRi using dCas9-Mxi.

## Guide-to-target assignments

The [`library_v1/sequence-good-targets.txt`](https://github.com/ingolia-lab/yeast-crispri/blob/master/library_v1/sequence-good-targets.txt) file provides a table
of the 61,094 guides in the Ingolia lab v1 guide library, along with
information about target genes. In this table,

* *Guide* gives a systematic name for each guide in the library, based
   on the systematic name of the gene target and a unique number. [See
   below](#degenerate-guides) for more information about degenerate guides marked
   `DEGEN`.

* *TargetLoc* gives the genomic location of the target sequence. This
   is listed in the format below, which is very similar to the
   information in BED-format annotations:

   _chromosome_:_start-end_(_strand_)

   where:
   * _chromosome_ is the reference sequence name   
   * _start_ is the starting position -- the 5'-most nucleotide in the
     target in 0-based coordinates.
   * _end_ is the ending position -- the 5'-most nucleotide *not* in
     the target in 0-based coordinates.
   * _strand_ is `+` or `-`, and a `+`-strand target will end in
     ...GG-3' while a `-`-strand target will start with 5'-CC.

* *Yorf1* is the systematic name of the primary, most likely
   target. This may be `NA` when a guide has only dubious gene
   targets, which are not candidates for *Yorf1* (or *Yorf2*,
   *Yorf3*).

* *Offset1* is the distance in nucleotides from the center of the
   guide sequence to the beginning of the primary target.

* *StartType1* indicates whether the primary target location is
   defined according to a transcription start site (`TSS`) or just the
   coding sequence (`CDS`).

* *Yorf2*, *Offset2*, and *StartType2* give the same information for a
   secondary, less likely target if one exists.

* *Yorf3*, *Offset3*, and *StartType3* give the same information for a
   tertiary, even less likely target if one exists (this is very
   rare).

* *Yorfs* gives an ordered, comma-separated list of targets

* *PrimaryGuide* gives the first choice for the guide systematic name,
   which generally matches *Guide* except in unusual circumstances
   involving dubious genes.

* *PrimaryGuideSource* indicates which of the potential targets was
   used to define the primary guide identity.

* *Oligo* gives the sequence of the guide oligo, with upstream and
   downstream constant sequence in lower-case and variable sequence in
   upper-case.

## Degenerate guides

A small fraction of guides are _degenerate_, meaning that they target
non-unique sequences that occur in multiple places in the
genome. These degenerate guides have `DEGEN` in their *Guide*
name. For these degenerate guides, the `sequence-good-targets.txt`
table lists information about one genomic target location. The
[`library_v1/sequence-redundant-targets.txt`](https://github.com/ingolia-lab/yeast-crispri/blob/master/library_v1/sequence-redundant-targets.txt) table gives similar
information about other genomic targets. The redundant targets table
has the same information as the main table of good targets, and with
an additional column *DegenGuide* that gives the guide name in the
main `sequence-good-targets.txt` table.