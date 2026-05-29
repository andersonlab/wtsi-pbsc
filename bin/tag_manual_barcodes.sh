#!/usr/bin/env bash
# Set rc:i:1 on reads whose CB tag matches a manual barcode list.
# Replaces an existing rc:i:0 in place; adds rc:i:1 if no rc tag is present.
# IMPORTANT: Reads that already carry rc:i:1 that are not in the list will be set to rc:i:0, overriding the existing tag.
#Barcodes file is a list of barcodes with no header


# Usage: tag_manual_barcodes.sh <barcodes_file> <in_bam> <out_bam>

set -euo pipefail

barcodes_file=$1
in_bam=$2
out_bam=$3

samtools view -h "$in_bam" \
| awk -v bf="$barcodes_file" '
    BEGIN {
        FS = OFS = "\t"
        while ((getline line < bf) > 0) {
            gsub(/[[:space:]]/, "", line)
            if (line != "") barcodes[line] = 1
        }
    }
    /^@/ { print; next }
    {
        cb = ""; rc_idx = 0; has_rc1 = 0
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^CB:Z:/) { split($i, a, ":"); cb = a[3] }
            else if ($i ~ /^rc:i:/) {
                rc_idx = i
                if ($i == "rc:i:1") has_rc1 = 1
            }
        }
        if (cb in barcodes) {
            if (rc_idx > 0) $rc_idx = "rc:i:1"   # replace existing rc:i:0 (or other)
            else            $0 = $0 OFS "rc:i:1" # no rc tag present, append
        } else {
            if (rc_idx > 0) $rc_idx = "rc:i:0"   # replace existing rc:i:0 (or other)
            else            $0 = $0 OFS "rc:i:0" # no rc tag present, append

        }
        print
    }' \
| samtools view -bo "$out_bam"
