#!/usr/bin/env python3
"""
Match Vireo donors to barcoded scRNAseq samples by barcode overlap.

Inputs:
  1) barcode_map.tsv : two columns
       <scRNAseq sample ID> <TAB> <Cell Ranger barcodes.tsv.gz>
     - barcodes.tsv.gz: one barcode per line (no header)

  2) vireo_map.tsv   : two columns
       <Vireo sample ID> <TAB> <Vireo donor_ids.tsv>
     - donor_ids.tsv: TSV with header; only first 2 columns are used:
         cell, donor_id
     - rows where donor_id is "doublet" or "unassigned" are ignored

Output (TSV to stdout):
  pool
  donor
  barcode_sample_id
  n_barcodes_in_common
  n_barcodes
  %_in_common
  rank

Notes:
  - Barcodes read from the Vireo files are reverse-complemented before
    matching.

Usage:
  python match-barcodes.py barcode_map.tsv vireo_map.tsv > match-barcodes.txt
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


_NON_AGCT = re.compile(r"[^AGCT]")
_RC_TABLE = str.maketrans({"A": "T", "C": "G", "G": "C", "T": "A"})


def clean_barcode(bc: str) -> str:
    """Remove any non-DNA character."""
    return _NON_AGCT.sub("", bc)


def revcomp_dna(bc: str) -> str:
    """Reverse complement."""
    return bc.translate(_RC_TABLE)[::-1]


def open_text_maybe_gz(path: str) -> io.TextIOBase:
    """Open a text file, transparently handling .gz compression."""
    p = Path(path)
    if p.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", newline="")
    return open(p, "r", encoding="utf-8", newline="")


def read_two_col_map(tsv_path: str) -> List[Tuple[str, str]]:
    """Read a 2-column TSV (no header assumed)."""
    pairs: List[Tuple[str, str]] = []
    with open(tsv_path, "r", encoding="utf-8", newline="") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"{tsv_path}:{line_no}: expected >=2 tab-separated columns")
            pairs.append((parts[0], parts[1]))
    return pairs


def load_barcode_index(barcode_map_tsv: str) -> Dict[str, List[str]]:
    """
    Build dict:
        barcode -> list of sample IDs (from the first input file)
    """
    barcode_to_samples: Dict[str, List[str]] = defaultdict(list)

    for sample_id, barcode_file in read_two_col_map(barcode_map_tsv):
        seen = set()
        with open_text_maybe_gz(barcode_file) as bf:
            for line in bf:
                bc = clean_barcode(line.strip())
                if not bc:
                    continue
                if bc in seen:
                    continue
                seen.add(bc)
                barcode_to_samples[bc].append(sample_id)

    return barcode_to_samples


def read_vireo_assignments(vireo_file: str) -> Dict[str, List[str]]:
    """
    Read a vireo assignment file into:
        donor_id -> list of barcodes (cell IDs)

      - Ignore header
      - Use only first two columns: cell, donor_id
      - Skip donor_id == "doublet" or "unassigned"
      - Reverse complement the barcodes
    """
    donor_to_barcodes: Dict[str, List[str]] = defaultdict(list)

    with open_text_maybe_gz(vireo_file) as f:
        first = True
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"{vireo_file}:{line_no}: expected >=2 tab-separated columns")

            if first:
                first = False
                if parts[0] == "cell" and parts[1] == "donor_id":
                    continue

            cell_clean = clean_barcode(parts[0])
            donor = parts[1]

            if donor in {"doublet", "unassigned"}:
                continue
            if not cell_clean:
                continue

            cell_rc = revcomp_dna(cell_clean)
            donor_to_barcodes[donor].append(cell_rc)

    return donor_to_barcodes


def dense_ranks_from_counts(counts: Dict[str, int]) -> Dict[str, int]:
    """
    Dense ranking by count.
    e.g. counts 10,10,9,7 -> ranks 1,1,2,3
    """
    unique_sorted = sorted(set(counts.values()), reverse=True)
    rank_by_count = {c: i + 1 for i, c in enumerate(unique_sorted)}
    return {k: rank_by_count[v] for k, v in counts.items()}


def main(argv: List[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Match barcodes by overlap.")
    ap.add_argument("barcode_map_tsv", help="2-col TSV: <sample_id> <barcodes.tsv.gz>")
    ap.add_argument("vireo_map_tsv", help="2-col TSV: <sample_id> <donor_ids.tsv>")
    args = ap.parse_args(argv)

    barcode_to_samples = load_barcode_index(args.barcode_map_tsv)

    out = sys.stdout
    out.write(
        "\t".join(
            [
                "pool",
                "donor",
                "barcode_sample_id",
                "n_barcodes_in_common",
                "n_barcodes",
                "%_in_common",
                "rank",
            ]
        )
        + "\n"
    )

    for pool, vireo_file in read_two_col_map(args.vireo_map_tsv):
        donor_to_barcodes = read_vireo_assignments(vireo_file)

        for donor, barcodes in donor_to_barcodes.items():
            unique_barcodes = set(barcodes)
            n_barcodes = len(unique_barcodes)
            if n_barcodes == 0:
                continue

            overlap_counts: Counter[str] = Counter()

            for bc in unique_barcodes:
                for sample_id in barcode_to_samples.get(bc, []):
                    overlap_counts[sample_id] += 1

            if not overlap_counts:
                continue

            ranks = dense_ranks_from_counts(dict(overlap_counts))

            for barcode_sample_id, n_common in sorted(
                overlap_counts.items(), key=lambda kv: (-kv[1], kv[0])
            ):
                pct_in_common = round((n_common / n_barcodes) * 100.0, 1)
                out.write(
                    "\t".join(
                        [
                            pool,
                            donor,
                            barcode_sample_id,
                            str(n_common),
                            str(n_barcodes),
                            f"{pct_in_common:.1f}",
                            str(ranks[barcode_sample_id]),
                        ]
                    )
                    + "\n"
                )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())