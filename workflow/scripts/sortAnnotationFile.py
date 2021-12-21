import pandas as pd
import numpy  as np
import click

# TODO add default value for feature
# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Variant/region output of lifover')
@click.option('--feature',
              'features',
              required=True,
              type=(str, float),
              multiple=True,
              help='All feature names we expect with additional default value')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file in tsv format.')
def cli(input_file, output_file, features):

    df = pd.read_csv(input_file, sep='\t')
    
    featureSet = set()
    for feature, default_value in features:
        featureSet.add(feature)
        if feature not in df.columns:
            df[feature] = default_value
        else:
            df[feature] = df[feature].fillna(value=default_value)
    cols = list(df.columns)
    for feature in list(df.columns):
        if feature not in featureSet:
            cols.remove(feature)
    out = df.fillna(0.0).reindex(["CHR", "POSITION", "ID"] + sorted(cols), axis=1)
    out.to_csv(output_file, sep='\t', index=None, header=True)


if __name__ == '__main__':
    cli()
