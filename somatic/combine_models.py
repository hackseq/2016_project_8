#!/usr/bin/env python
#
import pandas


tbl = {
    '0': "/hackseq/team8somatic/count_table_0pct.csv",
    '25': "/hackseq/team8somatic/count_table_25pct.csv",
    '50': "/hackseq/team8somatic/count_table_50pct.csv",
    '75': "/hackseq/team8somatic/count_table_75pct.csv",
    '100': "/hackseq/team8somatic/count_table_100pct.csv",
}

def combine_csvs(csv_labels):

    tbls = { name: pandas.DataFrame(fn) for (name, fn) in csv_labels }

    df = pandas.DataFrame()
    for (name, sub_df) in tbls.items():

        for col in sub_df.columns:
            if not col in ('chrom', 'pos'):
                df[col + "." + name] = sub_df

    df.to_csv("combined.csv")


if __name__ == "__main__":
    combine_csvs(tbl)