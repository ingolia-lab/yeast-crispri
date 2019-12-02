#!/bin/bash

set -x
set -e

if [[ -z "$1" ]];
then
    echo "Usage: $0 <DATADIR>"
    exit 1
fi

export DATADIR="$1"

./download_data.sh
cargo build
./target/debug/yeast_grna "${DATADIR}"
