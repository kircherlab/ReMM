import pandas as pd
import click
import vcfpy


# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Variant/region output of lifover')
@click.option('--regions',
              'regions_file',
              required=True,
              type=str,
              help='BED file with regions')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file in tsv format.')
def cli(input_file, regions_file, output_file):

    reader = vcfpy.Reader.from_path(input_file)
    regions = pd.read_csv(regions_file, sep='\t', header=None,comment="#")
    features = []
    for index, row in regions.iterrows():
        for record in reader.fetch(row[0], row[1], row[2]):
            df = pd.DataFrame(record.INFO)
            df['CHR'] = record.CHROM
            df['POSITION'] = record.POS
            if record.ID:
                df['ID'] = record.ID
            else:
                df['ID'] = "."
            features += [df]

    if len(features) > 1:
        out = pd.concat(features)
    elif len(features) == 1:
        out = features[0]
    else:
        raise Exception("Cannot find record on interval")

    cols = list(out.columns)
    for remove in ["CHR", "POSITION", "ID"]:
        cols.remove(remove)
    out = out.reset_index(drop=True).reindex(["CHR", "POSITION", "ID"] + cols, axis=1)
    out.to_csv(output_file, sep='\t', index=None, header=True, na_rep="NaN")


if __name__ == '__main__':
    cli()
