import pandas as pd
import click


# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Variant/region output of lifover')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file in tsv format.')
def cli(input_file, output_file):

    df = pd.read_csv(input_file, sep='\t')

    cols = list(df.columns)

    for remove in ["CHR", "POSITION", "ID"]:
        cols.remove(remove)
    out = df.fillna(0.0).reindex(["CHR", "POSITION", "ID"] + sorted(cols), axis=1)
    out.to_csv(output_file, sep='\t', index=None, header=True)


if __name__ == '__main__':
    cli()
