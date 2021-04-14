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
    parser = argparse.ArgumentParser(prog='jackhmmer')
    parser = argparse.ArgumentParser(description='Give this program a directory loaded with extracted orfs, to run'
                                                 'jackhmmer')
    parser.add_argument('-orf_dir', '--OrfDir', help='Give this command the directory of database ORFs')
    parser.add_argument('-seed', '--Seed', help='Give this command the directory of seed ORFs to use')
    parser.add_argument('-genetic_code', '--GeneCode',
                        help='Give this command the genetic code of the input sequences. Defaults to 11')
    args = parser.parse_args()
    return args

def main(args):
    original_path = os.getcwd()
    os.mkdir('HMMer')
    os.mkdir('HMMer/FEMALE')
    os.mkdir('HMMer/MALE')
    os.chdir(args.OrfDir)
    for file in sorted(os.listdir()):
        if file == 'FEMALE':
            os.chdir(original_path + '/' + args.OrfDir + '/' + file)
            for mito_file in sorted(os.listdir()):
                os.system('jackhmmer -N 10 --tblout {} {} {}' .format(os.path.join(original_path,'HMMer/FEMALE/', mito_file[:-4] + '_forf.tbl'),mito_file,os.path.join(original_path,args.Seed,'forfs.faa')))
        elif file == 'MALE':
            os.chdir(original_path + '/' + args.OrfDir + '/' + file)
            for mito_file in sorted(os.listdir()):
                os.system('jackhmmer -N 10 --tblout {} {} {}'.format(
                    os.path.join(original_path, 'HMMer/MALE/', mito_file[:-4] + '_morf.tbl'), mito_file,
                                 os.path.join(original_path, args.Seed, 'morfs.faa')))
        os.chdir(original_path + '/' + args.OrfDir)

if __name__ == "__main__":
    args = parse()
    main(args)