#!/usr/bin/env python3

import argparse
import pandas as pd


def remove_matching_rows(table1: pd.DataFrame, table2: pd.DataFrame) -> pd.DataFrame:
    # ensure required columns exist
    required = {"pool", "donor"}
    if not required.issubset(table1.columns):
        raise ValueError("table1 must contain 'pool' and 'donor' columns")
    if not required.issubset(table2.columns):
        raise ValueError("table2 must contain 'pool' and 'donor' columns")

    # mark matches
    merged = table1.merge(
        table2[["pool", "donor"]],
        on=["pool", "donor"],
        how="left",
        indicator=True
    )

    # keep only rows not in table2
    result = merged[merged["_merge"] == "left_only"].drop(columns="_merge")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Remove rows from table1 if (pool, donor) exists in table2"
    )
    parser.add_argument("table1", help="Input table1 CSV")
    parser.add_argument("table2", help="Input table2 CSV")
    parser.add_argument("output", help="Output CSV")

    args = parser.parse_args()

    df1 = pd.read_csv(args.table1)
    df2 = pd.read_csv(args.table2)

    result = remove_matching_rows(df1, df2)

    result.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
