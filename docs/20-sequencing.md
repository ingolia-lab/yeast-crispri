---
layout: page
title: Deep Sequencing
---

## Barcode-to-guide assignment

Primers NI-956 and NI-1038 carry out first-round PCR to create a
paired-end library that provides guide-to-barcode assignment. The
product from this first-round PCR should be amplified using primers
that contain the full Illumina P7 flowcell primer sequence, such as
the NEBNext Multiplex Oligos for Illumina.

The Illumina TruSeq Read 1 primer site is immediately adjacent to the
barcode, and reads the barcode for both assignment and counting
libraries. Primer NI-956 adds the Illumina P5 flowcell primer sequence
along with the rest of the sequencing primer binding site.

```
NI-956 5'-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC

    …TCAGGCGCGCCACTTCACGCATGCNNNNNNNNNNNNNNNNNNNNNNNNNAGATCGGAAGAGCGTCGTGCTATAGTGAGTCGTATTACATGCTCAAGAGCTCGATCCG…
NI-956 reverse complement                            <-GATCGGAAGAGCGTCGTG
                                                                         TAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
                                                                                            <- Illumina P5
```

Primer NI-1038 binds in the _P(RPR1)_ polymerase, 8 base pairs
upstream from the start of the sgRNA. It adds the Illumina TruSeq Read
2 primer site, and the product is then suitable for further
amplification with primers containing the Illumina P7 flowcell primer,
along with an index if needed. Because the final 8 bases of the
_P(RPR1)_ promoter are amplified from the plasmid library, any
mutations in these positions will be retained in the final sequencing
library.

```
NI-1038 5'-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTcgaaactctgggagctgc
              TruSeq Read #2 primer

                     …cgcggctgggaacgaaactctgggagctgcgattggcaCATTTCATTACCCGCAGAGCgtttt…
                                  cgaaactctgggagctgc->
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

## Barcode counting

In vitro transcription generates an RNA product beginning with the `G`
nucleotide at the end of the T7 promoter.

```
templ …AGCTATCAGGCGCGCCACTTCACGCATGCNNNNNNNNNNNNNNNNNNNNNNNNNAGATCGGAAGAGCGTCGTGCTATAGTGAGTCGTATTACATG…
product rev compl       …TTCACGCATGCNNNNNNNNNNNNNNNNNNNNNNNNNAGATCGGAAGAGCGTCGTGC
```

Reverse transcription with NI-1032 generates a fairly short DNA
product with a partial Illumina TruSeq Read 2 and flowcell P7 adapter
sequence. The product contains the partial Illumina TruSeq Read 1
primer site encoded by the plasmid library as well.

```
GACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNGCATGCGTGAAGTGGCGCGCCTGATAGCTCGTTTAAACTGGGTACCGGCCGCATAGCGAACGTGTAGGGCAGCGTTTCC…
NI-1032 reverse complement                                                  <-CTGGGTACCGGCCGCATA
                                                                                                AGATCGGAAGAGCACACGTCTGAA…
```

PCR using indexed primers such as the NEBNext Multiplex Oligos for
Illumina will amplify a sequencing library where the 25-base barcode
starts immediately at the beginning of read 1.

```
i5         AAT…ACAC(i5)ACACTCTTTCCCTACACGACGCTCTTCCGATCT
RT product reverse complement        CACGACGCTCTTCCGATCTNN…NNGCA…ATAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
i7 reverse complement                                               AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC(i7)AT…TG
```