import sys, os, subprocess
import argparse
from pandas import read_excel
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
Entrez.email='joao.edu.t@hotmail.com'
api_key = '105bc09b7b34973aa305f31be73aecbc4109'

def parse():
    """Parsing arguments"""
    parser = argparse.ArgumentParser(prog='getorf')
    parser = argparse.ArgumentParser(description='Give this program a directory loaded with mitogenomes, to run'
                                                 'the EMBOSS tool getorf')
    parser.add_argument('-mito', '--MitoFile', help='Give this command the directory of mitogenomes')
    parser.add_argument('-genetic_code', '--GeneCode',
                        help='Give this command the genetic code of the input sequences. Defaults to 11')
    args = parser.parse_args()
    return args

def main(args):
    original_path = os.getcwd()
    os.mkdir('Predicted_ORFs')
    os.mkdir('Predicted_ORFs/FEMALE')
    os.mkdir('Predicted_ORFs/MALE')
    os.chdir(args.MitoFile)
    for file in sorted(os.listdir()):
        if file == 'FEMALE':
            os.chdir(original_path + '/' + args.MitoFile + '/' + file)
            for mito_file in sorted(os.listdir()):
                os.system('getorf -sequence {} -table {} -minsize 30 -circular Yes -find 1 -methionine No -outseq {}' .format(mito_file,args.GeneCode,os.path.join(original_path,'Predicted_ORFs/FEMALE/',mito_file[:-3] + 'faa')))
        elif file == 'MALE':
            os.chdir(original_path + '/' + args.MitoFile + '/' + file)
            for mito_file in sorted(os.listdir()):
                os.system('getorf -sequence {} -table {} -minsize 30 -circular Yes -find 1 -methionine No -outseq {}'.format(
                        mito_file, str(args.GeneCode),
                        os.path.join(original_path, 'Predicted_ORFs/MALE/' + mito_file[:-3] + 'faa')))
        os.chdir(original_path + '/' + args.MitoFile)

if __name__ == "__main__":
    args = parse()
    main(args)