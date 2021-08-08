#!/bin/bash
set -euo pipefail

fasta_fn="$1"

if [[ -z $fasta_fn ]]; then
    >&2 echo "Usage: $0 /path/to/ref_genome.fa"
    exit 1
fi

dn="${fasta_fn}.unpack"
if [ "$#" -gt 1 ]; then
    dn="$2"
fi
if [[ -e $dn ]]; then
    >&2 echo "Destination already exists; to recreate, delete then retry: rm -rf '$dn'"
fi

mkdir -p "$dn"
dn=$(realpath "$dn")
seqkit fx2tab -i "$fasta_fn" | python3 -c "
import sys
import os
for line in sys.stdin:
    line = line.rstrip()
    sep = line.index('\t')
    fn = os.path.join('$dn',line[:sep])
    with open(fn, 'w') as outfile:
        outfile.write(line[sep+1:])
    print(fn)
    sys.stdout.flush()
"
echo "${dn}/"
