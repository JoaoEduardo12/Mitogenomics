import sys, os, subprocess
import argparse
from pandas import read_excel
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio.SeqRecord import SeqRecord

def parse():
    """Parsing arguments"""
    parser = argparse.ArgumentParser(prog='translate_file')
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-dir', '--Dir', help='Give this command the directory of database ORFs')
    parser.add_argument('-genetic_code', '--GeneCode',
                        help='Give this command the genetic code of the input sequences. Defaults to 11')
    args = parser.parse_args()
    return args

def main(args):
    os.chdir(args.Dir)
    for file in sorted(os.listdir()):
        with open(file[:-3] + 'faa','a+') as aa_fa:
            for dna_rec in SeqIO.parse(file,'fasta'):
                aa_fa.write('>{}\n' .format(dna_rec.id))
                aa_fa.write(str(dna_rec.seq.translate(table = args.GeneCode)) + '\n')
        aa_fa.close()

if __name__ == "__main__":
    args = parse()
    main(args)