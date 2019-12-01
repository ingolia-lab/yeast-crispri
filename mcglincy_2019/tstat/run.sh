#!/bin/bash

set -x
set -e

DATADIR="./data"
FIGUREDIR="./figures"

mkdir -p "${FIGUREDIR}"

export DATADIR
export FIGUREDIR
R --no-save < growth_tstat.R

