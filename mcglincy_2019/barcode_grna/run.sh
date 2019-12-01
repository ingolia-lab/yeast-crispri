#!/bin/bash

set -x
set -e

GUIDEDIR="../../library_v1/"
BARCODINGDIR="../../barcode_v1_181005/"

FIGUREDIR="./figures"

mkdir -p "${FIGUREDIR}"

export GUIDEDIR
export BARCODINGDIR
export FIGUREDIR
R --no-save < barcoding_stats.R
R --no-save < guide_stats.R

