import numpy as np

import pandas as pd
import click
import vcfpy


# options
@click.command()
@click.option('--folds',
              'folds_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Folds file which position belongs to which fold')
@click.option('--positives',
              'positives_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Annotated variant file with positives')
@click.option('--negatives',
              'negatives_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Annotated variant file with negatives')
@click.option('--output-data',
              'output_data_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file in tsv format.')
@click.option('--output-folds',
              'output_folds_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file in tsv format.')
@click.option('--output-labels',
              'output_label_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file in tsv format.')
def cli(folds_file, positives_file, negatives_file, output_data_file, output_folds_file, output_label_file):

    def loadData(f, label):
        df = pd.read_table(f, sep="\t")
        df["LABEL"] = label
        return df

    df = pd.concat([loadData(positives_file, 1), loadData(negatives_file, 0)])
    cytoband = pd.read_table(folds_file, header=None)

    df['FOLD'] = None
    cytoband[[1, 2]] = cytoband[[1, 2]]+1
    for i in range(len(cytoband)):
        ind = df[(df["CHR"] == cytoband[0][i]) & (cytoband[1][i] <= df["POSITION"]) & (cytoband[2][i] > df["POSITION"])].index
        if len(ind > 0):
            # print(len(ind))
            df.loc[ind, 'FOLD'] = cytoband[4][i]

    df['LABEL'].to_csv(output_label_file, sep='\t', header=None, index=None)
    df['FOLD'].to_csv(output_folds_file, sep='\t', header=None, index=None)

    cols = list(df.columns)
    for remove in ["CHR", "POSITION", "ID", "FOLD", "LABEL"]:
        cols.remove(remove)
    df[cols].to_csv(output_data_file, sep='\t', header=None, index=None, na_rep=0)


if __name__ == '__main__':
    cli()
