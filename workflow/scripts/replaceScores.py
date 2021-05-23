import click
import sys
import gzip


# options
@click.command()
# @click.option('--genome-score-file',
#               'genome_file',
#               required=True,
#               type=click.Path(exists=True, readable=True),
#               help='Whole genome tabix file with scores')
@click.option('--replace-score-file',
              'replace_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Scores to replace in the whole genome file. Must be sorted the same way!')
@click.option('--comment',
              'comments',
              required=False,
              type=str,
              multiple=True,
              help='Comment added on top of the file starting with #')
# def cli(genome_file, replace_file, comments):
def cli(replace_file, comments):

    # comments
    for comment in comments:
        print("#%s" % comment)
    # header
    print("#CHR\tPOS\tPROBABILITY")


    def getNextReplaceLine(replace_lines, i=-1):
        replace_end = False
        replace_line = None
        replace_split = None
        if (i+1 < len(replace_lines)):
            i += 1
            replace_line = replace_lines[i].strip()
            replace_split = replace_line.split("\t")
        else:
            replace_end = True
        return i, replace_line, replace_split, replace_end

    with gzip.open(replace_file, 'rt') as replace_f:
        replace_lines = replace_f.readlines()

        i, replace_line, replace_split, replace_end = getNextReplaceLine(replace_lines)

        replace_end = False 
        for genome_line in sys.stdin:
            genome_split = genome_line.strip().split("\t")
            contig = genome_split[0]
            position = int(genome_split[1])
            
            if not replace_end:

                # iterate until genome line is or can be equal
                while (contig == replace_split[0] and position > int(replace_split[1])):
                    print(replace_line)
                    i, replace_line, replace_split, replace_end = getNextReplaceLine(replace_lines, i)
                
                # print out replace line if == oherwise genome_line
                if contig == replace_split[0] and position == int(replace_split[1]):
                    print(replace_line)
                    i, replace_line, replace_split, replace_end = getNextReplaceLine(replace_lines, i)
                else:
                    print(genome_line.strip())
            else:
                print(genome_line.strip())


if __name__ == '__main__':
    cli()
