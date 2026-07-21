#!/bin/bash
# Extracts the unique IDs of a given feature type ("transcript" or "gene") from one or more GTF files.
# Usage: extract_gtf_feature_ids.sh <transcript|gene> <gtf_file> [<gtf_file> ...]
set -euo pipefail

feature_type="$1"
shift

awk -v ft="${feature_type}" -F'\t' '
  $3==ft {
    n_lines++
    attr=ft"_id \""
    if (match($9, attr"[^\"]+\"")) {
      s=substr($9,RSTART,RLENGTH)
      sub(attr,"",s)
      sub("\"$","",s)
      print s
      n_extracted++
    }
  }
  END {
    if ((n_lines+0) != (n_extracted+0)) {
      printf "ERROR: found %d %s feature line(s) but only extracted %d %s_id value(s) -- check for unquoted or malformed attributes\n", n_lines+0, ft, n_extracted+0, ft > "/dev/stderr"
      exit 1
    }
  }
' "$@" | sort -u
