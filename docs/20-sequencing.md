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

Primer NI-1038 binds in the _P(RPR1)_ promoter, 8 base pairs
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

## IVT-RT library generation protocol

* *Reagents*
    * Custom oligonucleotide NI-1032: `5'-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTATGCGGCCGGTACCCAG`
    * Zymoprep Yeast Plasmid Miniprep II (Zymo Research D2004)
    * XhoI (NEB R0146)
    * DNA Clean &amp; Concentrator-5 (Zymo Research D4013)
    * HiScribe T7 Quick High Yield RNA Synthesis Kit (NEB E2050S)
    * RNA Clean &amp; Concentrator-5 (Zymo Research R1013)
    * ProtoScript II (NEB M0368S)
    * Q5 High-Fidelity 2&times; Master Mix (NEB M0492S)
    * NEBNext Multiplex Oligos for Illumina (Dual Index Primer Set 1, NEB E7600S, or Dual Index Primer Set 2, NEB E7780S)
    * AMpure XP (Beckman Coulter A63880)
    * Agilent High Sensitivity DNA ScreenTape, reagents, and ladder

* *Equipment*
    * Magnetic separator for microfuge tubes
    * Agilent TapeStation

1. Harvest 1 - 5 &times; 10<sup>7</sup> yeast.

1. Resuspend in 1 ml sterile deionized water, pellet by centrifugation for 30 s at 10,000 &times; _g_, and aspirate all water.

   _Note:_ Washed yeast pellets can be stored at -80 ºC.

1. Extract DNA from yeast pellet using the Zymoprep Yeast Plasmid
   Miniprep II (Zymo Research D2004) according to the manufacturer's
   instructions, as described below.

    1. Add 200  µl "Solution 1" to the pellet

    1. Add 3.0  µl Zymolyase to resuspended pellet and mix by mild vortexing

    1. Incubate 1 hour at 37 ºC

    1. Add 200  µl "Solution 2" to the sample

    1. Add 400  µl "Solution 3" to the sample

    1. Centrifuge 3 min at maximum speed (16,000 - 20,000 &times; _g_)

    1. Transfer the supernatant to the Zymo-Spin-I column

    1. Spin the column 30 s at 10,000 - 16,000 &times; _g_ and discard the flow-through

    1. Add 550  µl "Wash Buffer", spin the column 1 min at 10,000 - 16,000 &times; _g_, and discard the flow-through

    1. Transfer the column to a new, clean microcentrifuge tube

    1. Add 10.0 µl Tris&bull;Cl 10 mM, pH 8.0 and centrifuge 1 min at 10,000 - 16,000 &times; _g_ to collect eluate containing extracted plasmid DNA

1. Linearize extracted DNA using XhoI (NEB R0146). Combine

   Volume  | Reagent                   | Final
   ------  | ------------------------- | -----
   10.0 µl | extracted plasmid DNA     |
   11.5 µl | deionized water           |
   2.5 µl  | 10&times; CutSmart buffer | 1&times;
   1.0 µl  | XhoI 20 U / µl            | 20 U

   Incubate 1 hour at 37 ºC

1. Purify linearized DNA using a DNA Clean &amp; Concentrator-5 (Zymo
   Research D4013) according to the manufacturer's instructions, with
   a 20 µl elution volume.

1. Prepare an _in vitro_ transcription reaction using the protocol for
   "short transcripts" in the HiScribe T7 Quick High Yield RNA
   Synthesis Kit (NEB E2050S). Combine

   Volume  | Reagent
   ------  | -------
   18.0 µl | purified, linearized DNA
   10.0 µl | NTP buffer mix
   2.0  µl | T7 RNA polymerase mix

   Incubate overnight at 37 ºC

1. Remove template DNA by adding 20.0 µl deionized water and 2.0 µl
   DNase I, 2 U / µl (provided in the kit) and incubate for 15 minutes at
   37 ºC.

1. Purify IVT product using an RNA Clean &amp; Concentrator-5 (Zymo Research R1013) according to the manufacturer's instructions, with a 15.0 µl elution volume.

1. Carry out reverse transcription using 10 ng of purified IVT product and ProtoScript II (NEB M0368S). Combine

   Quantity | Reagent
   -------- | -------
   10 ng    | purified IVT product RNA
   2.0  µl  | NI-1032 at 1 µM
   1.0  µl  | 10 mM ea. dNTPs
   _x_  µl  | water, to 10.0 µl final volume

   Denature 5 min at 65 ºC and place immediate on ice. Add

   Quantity | Reagent
   -------- | -------
   4.0  µl  | 5&times; ProtoScript II buffer
   2.0  µl  | 100 mM DTT
   1.0  µl  | ProtoScript II, 200 U / µl
   0.2  µl  | RNase inhibitor, 40 U / µl

   Incubate 1 hour at 42 ºC and then heat inactivate 20 min at 65 ºC.

1. Perform PCR amplification using Q5 master mix (NEB M0492S) and NEBNext Multiplex Oligos for Illumina (Dual Index Primer Set 1, NEB E7600S, or Dual Index Primer Set 2, NEB E7780S). Use distinct i5nn _and_ i7nn primers for each sample -- e.g., i501 and i701 for the first sample, i502 and i702 for the second sample, and so forth -- in order to exclude index hopping. Combine

   Quantity | Reagent
   -------- | -------
   5.0  µl  | Reverse transcription reaction
   2.5  µl  | NEBNext i5nn primer, 10 µM
   2.5  µl  | NEBNext i7nn primer, 10 µM
   15.0 µl  | deionized water
   25.0 µl  | Q5 High-Fidelity 2&times; Master Mix

   | Step                  |    | Temp  | Time
   | ----                  | -- | ----  | ----
   | Initial denaturation  |    | 98 ºC | 30 s
   | Amplification         | 7x | 98 ºC | 5 s
   |                       |    | 65 ºC | 10 s
   |                       |    | 72 ºC | 10 s
   | Final extension       |    | 72 ºC | 2 min

1. Purify PCR products using AMpure XP (Beckman Coulter A63880) according to the manufacturer's instructions, using 100 µl bead suspension for 50 µl PCR (a 2.0:1.0 bead-to-sample ratio). Elute DNA with 20. µl Tris&bull;Cl 10 mM, pH 8.0.

1. Quantify DNA by UV absorbance (e.g. on a Nanodrop spectrophotometer).

1. Validate the library using the High Sensitivity D1000 ScreenTape (Agilent), diluting as needed to ensure the ScreenTape sample is <1 ng / µl.

    * The final, dual-index library is 216 bp long, including 25 bp of barcode sequence and 55 bp of constant vector-derived sequence.

1. Carry out single-read, dual-indexed Illumina sequencing.

    * The first 25 bases of the sequencing read comprise the barcode, which has high nucleotide diversity. Constant plasmid-derived sequence follows the barcode, however, which is a monotemplate. 