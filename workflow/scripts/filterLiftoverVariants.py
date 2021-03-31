import pandas as pd
import click
from pyfaidx import Fasta
import difflib


# options
@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Variant/region output of lifover')
@click.option('--reference-old',
              'reference_old_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Reference of the old positions')
@click.option('--reference-new',
              'reference_new_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Reference of the new positions')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file in vcf format.')
def cli(input_file, reference_old_file, reference_new_file, output_file):

    reference_old = Fasta(reference_old_file,
                          one_based_attributes=False)
    reference_new = Fasta(reference_new_file,
                          one_based_attributes=False)

    var = pd.read_csv(input_file, sep='\t', header=None)

    ref_new = []
    diffs = []
    for i in range(len(var)):
        ref = reference_new[var[0][i]][var[1][i]:var[2][i]].seq
        n = 5
        ref_int_new = reference_new[var[0][i]][var[1][i]-n:var[2][i]+n].seq
        p = int(var[3][i].split(':')[1])
        ref_int_old = reference_old[var[0][i]][p-1-n:p+n].seq
        diff = sum([i[0] != ' ' for i in difflib.ndiff(ref_int_new, ref_int_old)]) / 2
        diffs.append(1-diff/1001)
        ref_new.append(ref)

    ref_old = var[3].str.split(':', expand=True)[2]
    vcf = pd.DataFrame()
    vcf[[0, 2]] = var[[0, 3]]
    vcf[1] = var[1]+1
    # ID

    # or remain
    # vcf[4]=var[3].str.split(':', expand = True)[3]
    vcf[3] = ref_new
    vcf[[5, 6]] = '.'
    vcf[4] = '.'

    vcf[7] = var[3].str.split(':', expand=True)[4]

    vcf = vcf[[0, 1, 2, 3, 4, 5, 6, 7]]

    vcf = vcf[vcf[2].str.split(':', expand=True)[3].str.len() == 1]
    vcf = vcf[vcf[2].str.split(':', expand=True)[2].str.len() == 1]

    vcf[3] = vcf[2].str.split(':', expand=True)[2]
    vcf[4] = vcf[2].str.split(':', expand=True)[3]

    vcf.columns = ["#CHROM", 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    out_file = open(output_file, 'a')
    out_file.write('##fileformat=VCFv4.1')

    vcf = vcf.T.reset_index().T

    vcf.to_csv(out_file, sep='\t', index=None)

    out_file.close()


if __name__ == '__main__':
    cli()
