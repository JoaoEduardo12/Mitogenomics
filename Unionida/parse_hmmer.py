import sys, os, subprocess
import argparse
from pandas import read_excel
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SearchIO

def parse():
    """Parsing arguments"""
    parser = argparse.ArgumentParser(prog='jackhmmer')
    parser = argparse.ArgumentParser(description='Give this program a directory loaded with extracted orfs, to run'
                                                 'jackhmmer')
    parser.add_argument('-orf_dir', '--OrfDir', help='Give this command the directory of database ORFs')
    parser.add_argument('-seq_dir', '--SeqDir', help='Give this command the directory of seed ORFs to use')
    parser.add_argument('-genetic_code', '--GeneCode',
                        help='Give this command the genetic code of the input sequences. Defaults to 11')
    args = parser.parse_args()
    return args

def main(args):
    original_path = os.getcwd()
    os.mkdir('ORF_hits')
    os.mkdir('ORF_hits/FEMALE')
    os.mkdir('ORF_hits/MALE')
    os.chdir(args.OrfDir)
    for file in sorted(os.listdir()):
        if file == 'FEMALE':
            os.chdir(original_path + '/' + args.OrfDir + '/' + file)
            for hmmFile in sorted(os.listdir()):
                with open(hmmFile,'rU') as handle:
                    for record in SearchIO.parse(handle, 'hmmer3-tab'):
                        query_id = record.id # seqID from fasta
                        hits = record.hits
                        num_hits = len(hits) #how many hits per query?
                        if num_hits > 0:
                            with open(os.path.join(original_path + '/ORF_hits/FEMALE', hmmFile[:-4] + '.faa'), 'a+') as final_file:
                                for dna_rec in SeqIO.parse(os.path.join(original_path, args.SeqDir + '/FEMALE', hmmFile[:-9] + '.faa'),'fasta'):
                                    if dna_rec.id == query_id:
                                        final_file.write('>{}\n{}\n' .format(dna_rec.description, str(dna_rec.seq)))
        elif file == 'MALE':
            os.chdir(original_path + '/' + args.OrfDir + '/' + file)
            for hmmFile in sorted(os.listdir()):
                with open(hmmFile,'rU') as handle:
                    for record in SearchIO.parse(handle, 'hmmer3-tab'):
                        query_id = record.id # seqID from fasta
                        hits = record.hits
                        num_hits = len(hits) #how many hits per query?
                        if num_hits > 0:
                            with open(os.path.join(original_path + '/ORF_hits/MALE', hmmFile[:-4] + '.faa'), 'a+') as final_file:
                                for dna_rec in SeqIO.parse(os.path.join(original_path, args.SeqDir + '/MALE', hmmFile[:-9] + '.faa'),'fasta'):
                                    if dna_rec.id == query_id:
                                        final_file.write('>{}\n{}\n' .format(dna_rec.description, str(dna_rec.seq)))
                            #for i in range(0,num_hits):
                             #   hmm_name = hits[i].id # hit name
                              #  hmm_description = hits[i].description # hit description
                               # current_evalue = hits[i].evalue #evalue of hit
                handle.close()



if __name__ == "__main__":
    args = parse()
    main(args)