import sys, os, subprocess
import argparse
from pandas import read_excel
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SearchIO
Entrez.email='joao.edu.t@hotmail.com'
api_key = '105bc09b7b34973aa305f31be73aecbc4109'

def parse():
    """Parsing arguments"""
    parser = argparse.ArgumentParser(prog='jackhmmer')
    parser = argparse.ArgumentParser(description='Give this program a directory loaded with extracted orfs, to run'
                                                 'jackhmmer')
    parser.add_argument('-dir', '--Dir', help='Give this command the directory of database ORFs')
    parser.add_argument('-genetic_code', '--GeneCode',
                        help='Give this command the genetic code of the input sequences. Defaults to 11')
    args = parser.parse_args()
    return args

def main(args):
    original_path = os.getcwd()
    os.mkdir('ORF_Alignments')
    os.mkdir('ORF_Alignments/FEMALE')
    os.mkdir('ORF_Alignments/MALE')
    os.chdir(args.Dir)
    for file in sorted(os.listdir()):
        if file == 'FEMALE':
            os.chdir(original_path + '/' + args.Dir + '/' + file)
            for seq_file in sorted(os.listdir()):
                with open(os.path.join(original_path,'ORF_Alignments/FEMALE/all_forfs.fas'), 'a+') as final_file:
                    with open(seq_file) as handle:
                        record = list(SeqIO.parse(handle,'fasta'))
                        final_file.write('>{}\n{}\n' .format(record[0].description,str(record[0].seq)))
                    handle.close()
                final_file.close()
        elif file == 'MALE':
            os.chdir(original_path + '/' + args.Dir + '/' + file)
            for seq_file in sorted(os.listdir()):
                with open(os.path.join(original_path,'ORF_Alignments/MALE/all_morfs.fas'), 'a+') as final_file:
                    with open(seq_file) as handle:
                        record = list(SeqIO.parse(handle,'fasta'))
                        final_file.write('>{}\n{}\n' .format(record[0].description,str(record[0].seq)))
                    handle.close()
                final_file.close()
    os.chdir(original_path + '/ORF_Alignments/FEMALE')
    os.system('t_coffee -seq all_forfs.fas')
    os.chdir(original_path + '/ORF_Alignments/MALE')
    os.system('t_coffee -seq all_morfs.fas')
    os.system('hmmsearch --max --tblout output.hmm F_profile.hmm /home/edu/Desktop/Bioinformatica/Unionida/NCBI_ORFs/forfs.faa')

if __name__ == "__main__":
    args = parse()
    main(args)