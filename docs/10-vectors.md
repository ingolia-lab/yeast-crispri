---
layout: page
title: Vectors
---

## CRISPRi effector protein
*pNTI647* Expression of the CRISPRi effector protein dCas9-Mxi1, along
with the TetR regulator controlling guide expression, is provided by
pNTI647. This integrating plasmid is based on the [KanR-marked
EasyClone 2.0 vector pCfB2225](https://www.addgene.org/67553/).

## Inducible guide RNAs
*pNTI661* The base guide RNA expression vector is an episomal,
Leu-marked plasmid. Barcoded guide RNA libraries are constructed in a
two-step process.

Guide RNAs are introduced by digesting the pNTI661 vector with BamHI
and HindIII, amplifying the guide RNA pool using primers NM636 and
NM637, and performing Gibson-style assembly (we like NEBuilder HiFi)
to introduce guide RNAs into the plasmid:

```
pNTI661     …cgcggctgggaacgaaactctgggagctgcgattgg/gatcctcgaagct/ttttagagctagaaatagcaagttaaaataaggctag…

digest  …cgcggctgggaacgaaactctgggagctgcgattgg     gatcctcgaagct     ttttagagctagaaatagcaagttaaaataaggctag…
            |||||||||||||||||||||||||||||||||                       ||||||||||||||||||||||||||||||||||
PCR      5'-ggctgggaacgaaactctgggagctgcgattggcaCATTTCATTACCCGCAGAGCgttttagagctagaaatagcaagttaaaataaggc-3'

NM636       ggctgggaacgaaactctgggagctgcgattggca
pool                       tctgggagctgcgattggcaCATTTCATTACCCGCAGAGCgttttagagctagaaatagc
NM637 reverse complement                                           gttttagagctagaaatagcaagttaaaataaggc
```

After the guide library is cloned and amplified, barcodes are
introduced in a second Gibson-style assembly reaction. The library is
digested with SphI, leaving 3' overhangs that are not resected during
the assembly reaction. The barcodes are amplified from a degenerate
oligo template using primers that introduce the Illumina TruSeq Read 1
primer site and a T7 RNA polymerase promoter.

```
pNTI661                                 …AGCTATCAGGCGCGCCACTTCACG/CATG/CTCAAGAGCTCGATCCGCAGGC…    

digest …AGCTATCAGGCGCGCCACTTCACGCATG                                                               CATGCTCAAGAGCTCGATCCGCAGGC…
           |||||||||||||||||||||||||                                                               |||||||||||||||||||||||
PCR     5'-TATCAGGCGCGCCACTTCACGCATGCNNNNNNNNNNNNNNNNNNNNNNNNNAGATCGGAAGAGCGTCGTGCTATAGTGAGTCGTATTACATGCTCAAGAGCTCGATCCGCA-3'

NI-1027    TATCAGGCGCGCCACTTCACGCATGC
NI-1026             CGCCACTTCACGCATGCNNNNNNNNNNNNNNNNNNNNNNNNNAGATCGGAAGAGCGTCGT
NI-1041 reverse complement                                    AGATCGGAAGAGCGTCGTGCTATAGTGAGTCGTATTACATGCTCAAGAGCTCGATCCGCA
```
