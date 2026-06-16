#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: append_to_tag.sh -t TARGET_TAG (-s STRING | -T TAG) [-i <in.bam|->] [-o <out.bam|->] [-O FMT] [-d DELIM] [-p CPUS] [-n]

Appends a string to the VALUE of a chosen SAM tag (e.g. CB:Z).
Reads BAM/SAM/CRAM from a file or stdin; writes SAM/BAM to a file or stdout.

The appended value is one of:
  -s STRING   a fixed string (e.g. a sample id), OR
  -T TAG      the value of another SAM tag from the same read

Options:
  -t TAG    target tag to modify                    (required), e.g. CB:Z
  -i FILE   input BAM, or - for stdin               (default: -)
  -o FILE   output, or - for stdout                 (default: -)
  -O FMT    output format: sam | bam | ubam          (default: bam)
            sam  = uncompressed text, ideal for piping into another of these scripts
            bam  = compressed BAM (indexable)
            ubam = uncompressed BAM
  -s STR    fixed string to append                  (mutually exclusive with -T)
  -T TAG    source tag whose value to append        (mutually exclusive with -s)
  -d CHAR   delimiter between tag value and value     (default: _)
  -p INT    threads for samtools                      (default: 1)
  -n        do not index the output                  (auto-skipped when writing to stdout)
  -h        show this help

Notes:
  * Reads lacking the target tag are left unchanged (the tag is NOT created).
  * With -T, if the source tag is absent the target is left unchanged.
  * Header lines (^@) are passed through untouched.
  * -O sam writes the awk stream straight to output with no BAM (de)coding;
    only bam/ubam output can be indexed.
  * Indexing needs a coordinate-sorted BAM/uBAM file; use -n for unsorted file output.
  * To drop supplementary alignments, add `-F 0x800` to the first `samtools view`.
EOF
  exit "${1:-1}"
}

tgt_tag="" in_bam="-" out_bam="-" out_fmt="bam" fixed="" src_tag="" delim="_" cpus=1 no_index=0

while getopts ":t:i:o:O:s:T:d:p:nh" opt; do
  case $opt in
    t)  tgt_tag=$OPTARG ;;
    i)  in_bam=$OPTARG ;;
    o)  out_bam=$OPTARG ;;
    O)  out_fmt=$OPTARG ;;
    s)  fixed=$OPTARG ;;
    T)  src_tag=$OPTARG ;;
    d)  delim=$OPTARG ;;
    p)  cpus=$OPTARG ;;
    n)  no_index=1 ;;
    h)  usage 0 ;;
    :)  echo "Error: -$OPTARG requires an argument" >&2; usage ;;
    \?) echo "Error: unknown option -$OPTARG" >&2; usage ;;
  esac
done

[[ -n $tgt_tag ]] || { echo "Error: -t (target tag) is required" >&2; usage; }

# require exactly one source
if { [[ -n $fixed ]] && [[ -n $src_tag ]]; } || { [[ -z $fixed ]] && [[ -z $src_tag ]]; }; then
  echo "Error: provide exactly one of -s (fixed string) or -T (source tag)" >&2
  usage
fi

mode="string"; src="$fixed"
[[ -n $src_tag ]] && { mode="tag"; src="$src_tag"; }

# transform() emits a valid SAM stream (header + tab-delimited records) on stdout
transform() {
  samtools view -@ "$cpus" -h "$in_bam" \
  | awk -v tgt="$tgt_tag" -v mode="$mode" -v src="$src" -v delim="$delim" 'BEGIN{FS=OFS="\t"}
  {
    if ($0 ~ /^@/) { print; next }

    # 1) resolve the value to append
    app = ""
    if (mode == "string") {
      app = src
    } else {
      sprefix = src ":"; splen = length(sprefix)
      for (i=12; i<=NF; i++) {
        if (substr($i, 1, splen) == sprefix) { app = substr($i, splen+1); break }
      }
    }

    # 2) append it to the target tag if both the value and the tag are present
    if (app != "") {
      tprefix = tgt ":"; tplen = length(tprefix)
      for (i=12; i<=NF; i++) {
        if (substr($i, 1, tplen) == tprefix) { $i = $i delim app; break }
      }
    }
    print
  }'
}

# route to the requested output format
case "$out_fmt" in
  sam)  # awk already produces SAM; no re-encode needed
    if [[ "$out_bam" == "-" ]]; then transform; else transform > "$out_bam"; fi ;;
  ubam) transform | samtools view -@ "$cpus" -h -u -o "$out_bam" - ;;
  bam)  transform | samtools view -@ "$cpus" -h -b -o "$out_bam" - ;;
  *)    echo "Error: -O must be one of sam | bam | ubam" >&2; usage ;;
esac

# index only a real (non-stdout) BAM/uBAM file
if [[ "$out_fmt" != "sam" && "$out_bam" != "-" && "$no_index" -eq 0 ]]; then
  samtools index -@ "$cpus" "$out_bam"
fi