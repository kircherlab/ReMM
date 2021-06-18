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
              type=str,
              multiple=True,
              help='All feature names we expect')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file in tsv format.')
def cli(input_file, output_file, features):

    df = pd.read_csv(input_file, sep='\t')

    

    for feature in features:
        if feature not in df.columns:
            df[feature] = np.nan
    cols = list(df.columns)
    for feature in list(df.columns):
        if feature not in features:
            cols.remove(feature)
    out = df.fillna(0.0).reindex(["CHR", "POSITION", "ID"] + sorted(cols), axis=1)
    out.to_csv(output_file, sep='\t', index=None, header=True)


if __name__ == '__main__':
    cli()
