#!/usr/bin/env python3

import argparse
import pandas as pd


def replace_from_table2(table1: pd.DataFrame, table2: pd.DataFrame) -> pd.DataFrame:
    # check columns match
    if set(table1.columns) != set(table2.columns):
        raise ValueError(
            f"Tables have different columns.\n"
            f"table1: {sorted(table1.columns)}\n"
            f"table2: {sorted(table2.columns)}"
        )

    # columns that may need replacing
    replace_cols = [c for c in ["donor_gt", "barcode_sample_id"] if c in table1.columns]

    # nothing to replace
    if not replace_cols:
        return table1.copy()

    # merge
    merged = table1.merge(
        table2[["pool", "donor"] + replace_cols],
        on=["pool", "donor"],
        how="left",
        suffixes=("", "_new")
    )

    # replace values
    for col in replace_cols:
        merged[col] = merged[f"{col}_new"].combine_first(merged[col])
        merged.drop(columns=f"{col}_new", inplace=True)

    return merged


def main():
    parser = argparse.ArgumentParser(description="Replace donor columns using pool/donor key")
    parser.add_argument("table1", help="First CSV file (to modify)")
    parser.add_argument("table2", help="Second CSV file (source values)")
    parser.add_argument("output", help="Output CSV file")

    args = parser.parse_args()

    # read CSVs
    df1 = pd.read_csv(args.table1)
    df2 = pd.read_csv(args.table2)

    # process
    result = replace_from_table2(df1, df2)

    # save
    result.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
