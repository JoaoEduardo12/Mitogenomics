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
    parser = argparse.ArgumentParser(prog='orf_analysis')
    parser = argparse.ArgumentParser(description='Give this program the search or database arguments'
                                                 ' to search and find orfs')
    parser.add_argument('-search', '--Search',
                        help='Give this command a file with the search term in NCBI '
                             '(right side of the page in "Search Details, when you query an entry")')
    parser.add_argument('-taxonomy', '--Tax',
                        help='Give this command the highest recent common taxonomy classification '
                             'within search (optional)')
    parser.add_argument('-mito', '--MitoFile', help='Give this command the directory of mitogenomes')
    parser.add_argument('-ann', '--AnnFile', help='Give this command the directory of annotations')
    parser.add_argument('-csv', '--SpSeq',
                        help='Give this command the Species Sequence excel file, to get taxonomy information'
                             ' (optional)')
    parser.add_argument('-genetic_code', '--GeneCode',
                        help='Give this command the genetic code of the input sequences. Defaults to 11')
    args = parser.parse_args()
    return args

def main(args):
    original_path = os.getcwd()
    os.mkdir('ORFs')
    os.chdir(args.MitoFile)
    for file in sorted(os.listdir()):
        if file == 'FEMALE' or file == 'MALE':
            if file == 'FEMALE':
                sex_type = 'F'
            elif file == 'MALE':
                sex_type = 'M'
            os.chdir(original_path + '/' + args.MitoFile + '/' + file)
            for mito_file in sorted(os.listdir()):
                mitogenome, annotation = read_mitogenomes(mito_file, original_path, file)
                extract_orfs(mitogenome,annotation, original_path)
    get_ncbi_orfs(args, original_path)

def read_mitogenomes(mito_file, path, file):
    print('Extracting info on {}' .format(mito_file[:-4]), end = '\r')
    mitogenome = SeqIO.read(mito_file,'fasta')
    os.chdir(path + '/' + args.AnnFile + '/' + file)
    for ann_file in sorted(os.listdir()):
        if ann_file.endswith('geneAnnotation.tab'):
            if ann_file[:-19] == mito_file[:-4]:
                annotation = open(ann_file, 'r')
        elif ann_file.endswith('geneAnnotations.tab'):
            if ann_file[:-20] == mito_file[:-4]:
                annotation = open(ann_file,'r')
    os.chdir(path + '/' + args.MitoFile + '/' + file)
    return mitogenome, annotation

def extract_orfs(mitogenome, annotation, path):
    mitogenome_rc = mitogenome.reverse_complement()
    mitogenome_sequence = str(mitogenome.seq)
    mitogenome_rc = str(mitogenome_rc.seq)
    mitogenome_sequence = mitogenome_sequence.upper()
    mitogenome_rc = mitogenome_rc.upper()
    for lines in annotation:
        if '---' not in lines and ' ' not in lines:
            molecule, type_, gene, start, stop, strain = lines.split('\t')
            if gene == 'FORF':
                if strain == '1\n':
                    file = open(os.path.join(path, 'ORFs/forfs.fas'), 'a+')
                    file.write('>{}\n{}\n' .format(molecule,mitogenome_sequence[int(start)-1:int(stop)]))
                    file.close()
                elif strain == '-1\n':
                    file = open(os.path.join(path, 'ORFs/forfs.fas'), 'a+')
                    file.write('>{}\n{}\n'.format(molecule, mitogenome_rc[int(start)-1:int(stop)]))
                    file.close()
            if gene == 'MORF':
                if strain == '1\n':
                    file = open(os.path.join(path, 'ORFs/morfs.fas'), 'a+')
                    file.write('>{}\n{}\n' .format(molecule,mitogenome_sequence[int(start)-1:int(stop)]))
                    file.close()
                elif strain == '-1\n':
                    file = open(os.path.join(path, 'ORFs/morfs.fas'), 'a+')
                    file.write('>{}\n{}\n'.format(molecule, mitogenome_rc[int(start)-1:int(stop)]))
                    file.close()
            if gene == 'HORF':
                if strain == '1\n':
                    file = open(os.path.join(path, 'ORFs/horfs.fas'), 'a+')
                    file.write('>{}\n{}\n' .format(molecule,mitogenome_sequence[int(start)-1:int(stop)]))
                    file.close()
                elif strain == '-1\n':
                    file = open(os.path.join(path, 'ORFs/horfs.fas'), 'a+')
                    file.write('>{}\n{}\n'.format(molecule, mitogenome_rc[int(start)-1:int(stop)]))
                    file.close()

def get_ncbi_orfs(args, path):
    os.chdir(path)
    os.mkdir('NCBI_ORFs')
    print('\nSearching in NCBI... Hold tight!\n')
    handle = Entrez.esearch(db='nucleotide', retmax = 100000, term = args.Search, idtype = 'acc')
    record = Entrez.read(handle)
    print('\nThere are {} results found\n' .format(str(len(record['IdList']))))
    handle.close()
    forfs = ['forf','Forf','FOrf','FORf','FORF','fOrf','fORf','fORF','FoRf','fOrF','f orf','F orf','F Orf','F ORf','F ORF','f Orf','f ORf','f ORF','F oRf','f OrF','f-orf','F-orf','F-Orf','F-ORf','F-ORF','f-Orf','f-ORf','f-ORF','F-oRf','f-OrF','f_orf','F_orf','F_Orf','F_ORf','F_ORF','f_Orf','f_ORf','f_ORF','F_oRf','f_OrF']
    morfs = ['morf', 'Morf', 'MOrf', 'MORf', 'MORF', 'mOrf', 'mORf', 'mORF', 'MoRf', 'mOrF', 'm orf', 'M orf', 'M Orf',
             'M ORf', 'M ORF', 'm Orf', 'm ORf', 'm ORF', 'M oRf', 'm OrF', 'm-orf', 'M-orf', 'M-Orf', 'M-ORf', 'M-ORF',
             'm-Orf', 'm-ORf', 'm-ORF', 'M-oRf', 'm-OrF', 'm_orf', 'M_orf', 'M_Orf', 'M_ORf', 'M_ORF', 'm_Orf', 'm_ORf',
             'm_ORF', 'M_oRf', 'm_OrF']
    horfs = ['horf', 'Horf', 'HOrf', 'HORf', 'HORF', 'hOrf', 'hORf', 'hORF', 'HoRf', 'hOrF', 'h orf', 'H orf', 'H Orf',
             'H ORf', 'H ORF', 'h Orf', 'h ORf', 'h ORF', 'H oRf', 'h OrF', 'h-orf', 'H-orf', 'H-Orf', 'H-ORf', 'H-ORF',
             'h-Orf', 'h-ORf', 'h-ORF', 'H-oRf', 'h-OrF', 'h_orf', 'H_orf', 'H_Orf', 'H_ORf', 'H_ORF', 'h_Orf', 'h_ORf',
             'h_ORF', 'H_oRf', 'h_OrF']
    print('\nExtracting ORFs!\n')
    i = 0
    for rec_name in record['IdList']:
        i += 1
        print(' ' + str(i) + '/' + str(len(record['IdList'])), end='\r')
        #mitogenome = Entrez.efetch(db = 'nucleotide', id = rec_name, rettype = 'fasta', retmode = 'text')
        #mitogenome_sequence = mitogenome.read()
        with Entrez.efetch(db = 'nucleotide', rettype = 'gb', retmode = 'text', id = rec_name) as handle:
            for rec in SeqIO.parse(handle,'gb'):
                if rec.features:
                    for feature in rec.features:
                        if feature.type == 'CDS':
                            if 'gene' in feature.qualifiers.keys():
                                if feature.qualifiers['gene'][0].strip('') in forfs:
                                    file = open(os.path.join(path, 'NCBI_ORFs/forfs.fas'), 'a+')
                                    file.write('>{}\n{}\n'.format(rec.id, feature.location.extract(rec).seq))
                                    file.close()
                                elif feature.qualifiers['gene'][0].strip('') in morfs:
                                    file = open(os.path.join(path, 'NCBI_ORFs/morfs.fas'), 'a+')
                                    file.write('>{}\n{}\n'.format(rec.id, feature.location.extract(rec).seq))
                                    file.close()
                                elif feature.qualifiers['gene'][0].strip('') in horfs:
                                    file = open(os.path.join(path, 'NCBI_ORFs/horfs.fas'), 'a+')
                                    file.write('>{}\n{}\n'.format(rec.id, feature.location.extract(rec).seq))
                                    file.close()
                            else:
                                if feature.qualifiers['product'][0].strip('') in forfs:
                                    file = open(os.path.join(path, 'NCBI_ORFs/forfs.fas'), 'a+')
                                    file.write('>{}\n{}\n'.format(rec.id, feature.location.extract(rec).seq))
                                    file.close()
                                elif feature.qualifiers['product'][0].strip('') in morfs:
                                    file = open(os.path.join(path, 'NCBI_ORFs/morfs.fas'), 'a+')
                                    file.write('>{}\n{}\n'.format(rec.id, feature.location.extract(rec).seq))
                                    file.close()
                                elif feature.qualifiers['product'][0].strip('') in horfs:
                                    file = open(os.path.join(path, 'NCBI_ORFs/horfs.fas'), 'a+')
                                    file.write('>{}\n{}\n'.format(rec.id, feature.location.extract(rec).seq))
                                    file.close()

if __name__ == "__main__":
    args = parse()
    main(args)