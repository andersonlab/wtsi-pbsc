#!/bin/bash
# Usage: cache_isoquant_refs.sh <genedb> <fasta>
# Copies the IsoQuant reference DB and genome FASTA to node-local /tmp once per
# node (flock-protected, atomic via tmp+mv). Prints DB_LOCAL and FA_LOCAL paths
# to stdout on separate lines so the caller can capture them.
set -euo pipefail

genedb="$1"
fasta="$2"

ISOQUANT_CACHE="/tmp/isoquant_ref_cache"
mkdir -p "${ISOQUANT_CACHE}"

DB_LOCAL="${ISOQUANT_CACHE}/$(basename "${genedb}")"
FA_LOCAL="${ISOQUANT_CACHE}/$(basename "${fasta}")"

(
  flock -x 200
  if [ ! -f "${DB_LOCAL}" ]; then
    cp "${genedb}" "${DB_LOCAL}.tmp" && mv "${DB_LOCAL}.tmp" "${DB_LOCAL}"
  fi
) 200>"${ISOQUANT_CACHE}/.lock_db"

(
  flock -x 201
  if [ ! -f "${FA_LOCAL}" ]; then
    cp "${fasta}" "${FA_LOCAL}.tmp" && mv "${FA_LOCAL}.tmp" "${FA_LOCAL}"
  fi
  if [ ! -f "${FA_LOCAL}.fai" ]; then
    cp "${fasta}.fai" "${FA_LOCAL}.fai.tmp" && mv "${FA_LOCAL}.fai.tmp" "${FA_LOCAL}.fai"
  fi
) 201>"${ISOQUANT_CACHE}/.lock_fa"

echo "${DB_LOCAL}"
echo "${FA_LOCAL}"
