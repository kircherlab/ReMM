import click
import vcfpy
import pandas as pd


# options
@click.command()
@click.option('--feature-set',
              'feature_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='VCF feature set file')
@click.option('--variants',
              'variant_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Variants to annotate in VCF format (only positions)')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file in tsv format.')
def cli(feature_file, variant_file, output_file):

    variant_reader = vcfpy.Reader.from_path(variant_file)
    variants = {}
    for variant_record in variant_reader:
        if variant_record.CHROM not in variants:
            variants[variant_record.CHROM] = {}
        if variant_record.ID:
            variants[variant_record.CHROM][variant_record.POS] = variant_record.ID
        else:
            variants[variant_record.CHROM][variant_record.POS] = "."

    feature_reader = vcfpy.Reader.from_path(feature_file)
    features = []
    chrs = set()
    for record in feature_reader:
        # just for logging
        if record.CHROM not in chrs:
            print(record.CHROM)
            chrs.add(record.CHROM)

        if record.CHROM in variants and record.POS in variants[record.CHROM]:
            df = pd.DataFrame(record.INFO)
            df['CHR'] = record.CHROM
            df['POSITION'] = record.POS
            df['ID'] = variants[variant_record.CHROM][variant_record.POS]
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
