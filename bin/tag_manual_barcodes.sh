#!/usr/bin/env bash
# Add rc:i:1 tag to reads whose CB tag matches a manual barcode list.
# Reads that already carry rc:i:1 are left untouched.
#
# Usage: tag_manual_barcodes.sh <barcodes_file> <in_bam> <out_bam>

set -euo pipefail

barcodes_file=$1
in_bam=$2
out_bam=$3

samtools view -h "$in_bam" \
| awk -v bf="$barcodes_file" '
    BEGIN {
        while ((getline line < bf) > 0) {
            gsub(/[[:space:]]/, "", line)
            if (line != "") barcodes[line] = 1
        }
    }
    /^@/ { print; next }
    {
        cb = ""; has_rc1 = 0
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^CB:Z:/) { split($i, a, ":"); cb = a[3] }
            if ($i == "rc:i:1") has_rc1 = 1
        }
        if (cb in barcodes && !has_rc1) print $0 "\trc:i:1"
        else print
    }' \
| samtools view -bo "$out_bam"
