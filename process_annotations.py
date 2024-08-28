import click as ck
import json
from pathlib import Path

@ck.command()
@ck.option('--annots-file', '-af', required=True, help='Gene annnotations file')
@ck.option('--output-file', '-of', required=True, help='Output file')
def main(annots_file, output_file):
    with open(annots_file) as f:
        annots = {}
        next(f)
        with open(annots_file) as f:
            for line in f:
                it = line.strip().split('\t')
                gene = it[0]
                hp_term = it[2]
                if gene not in annots:
                    annots[gene] = []
                annots[gene].append(hp_term)

    with open(output_file, 'w') as f:
        for gene, terms in annots.items():
            f.write(' '.join(terms) + '\n')

if __name__ == '__main__':
    main()
