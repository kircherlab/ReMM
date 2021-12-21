import numpy as np
from pyfaidx import Fasta
import pandas as pd
import click

# options


@click.command()
@click.option('--input',
              'input_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Reference genoem fasta file')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output bed file of not N positions.')
def cli(input_file, output_file):

    contigs = ["chr%s" % str(c) for c in list(range(1, 23)) + ["Y", "X"]]

    f = Fasta(input_file, rebuild=False)

    output = []

    for contig in contigs:
        start_position = -1
        position = -1
        defined = np.where(np.asarray(f[contig]) != b'N')[0]
        for i in defined:
            if start_position < 0:
                start_position = i
                position = i
            elif position + 1 == i:
                position = i
            else:
                output += [[contig, start_position, position+1]]
                start_position = i
                position = i
        output += [[contig, start_position, position+1]]

    pd.DataFrame(output).to_csv(output_file, sep='\t', index=None, header=None)


if __name__ == '__main__':
    cli()
