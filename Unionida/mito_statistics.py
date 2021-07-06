import sys, os, subprocess
import argparse
from pandas import read_excel
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio.SeqUtils import GC
from CAI import RSCU, CAI, relative_adaptiveness
from Bio.SeqUtils import CodonUsage as CU
from Bio import AlignIO
Entrez.email='joao.edu.t@hotmail.com'

def parse():
	parser = argparse.ArgumentParser(prog="gene_stats")
	parser = argparse.ArgumentParser(description='Give both the directories of where the mitogenomes are and where the annotation tables are')
	parser.add_argument('-mito', '--MitoFile', help='Give this command the directory of mitogenomes')
	parser.add_argument('-ann', '--AnnFile', help='Give this command the directory of annotations')
	parser.add_argument('-csv', '--SpSeq', help ='Give this command the Species Sequence excel file, to get taxonomy information (optional)')
	parser.add_argument('-genetic_code', '--GeneCode', help = 'Give this command the genetic code of the input sequences. Defaults to 11')
	parser.add_argument('-compare', '--CompGen', help = '(Optional) Give this command two sets of taxonomic ranks, family, subfamily, tribe, even species names, or sex (in freshwater mussels), example: Mulleridae x Unionidae; Ambleminae x Anodontinae; Male x Female. This is important to perform several evolutionary analyses on a specific set of sequences, otherwise it performs All vs All')
	parser.add_argument('-CAI', '--Cai', help = '(Optional) Give this command the samples you wish to calculate CAI index on PCGs, e.g., F, M, Family, Subfamily, Tribe')
	args = parser.parse_args()
	return args

def main(args):
	print('\n')
	original_path = os.getcwd()
	os.mkdir('Mitogenome_Statistics')
	os.chdir(args.MitoFile)
	total_gene_dict = {}
	total_mitogenome_dict = {}
	for file in sorted(os.listdir()):
		if file == 'FEMALE' or file == 'MALE':
			if file == 'FEMALE':
				sex_type = 'F'
			elif file == 'MALE':
				sex_type = 'M'
			os.chdir(original_path + '/' + args.MitoFile + '/' + file)
			for mito_file in sorted(os.listdir()):
				mitogenome, annotation = read_mitogenomes(mito_file, original_path, file)
				mito_length = genome_size(mitogenome)
				mito_gc = mitogenome_gc_content(mitogenome)
				p_g, p_c, p_a, p_t, gc_sk, at_sk, p_g3, p_c3, p_a3, p_t3, gc_sk3, at_sk3 = calculate_skew(mitogenome) #???
				gene_characteristics, current_go, gene_dictionary = gene_calculations(mitogenome,annotation,args,mito_file)
				total_gene_dict[mito_file[:-4]] = gene_dictionary
				total_mitogenome_dict[mito_file[:-4]] = mitogenome
				mitogenome_order = gene_orders(current_go)
				if args.SpSeq:
					family, subfamily, tribe, species_name = read_taxonid(args.SpSeq, mito_file, sex_type, original_path)
				write_output(file = mito_file, 
					content =  'Mitogenome', 
					family = family, 
					subfamily = subfamily, 
					tribe = tribe, 
					species_name = species_name, 
					length = mito_length, 
					gc = mito_gc, 
					sex = sex_type, 
					gene_char = gene_characteristics, 
					order = mitogenome_order, 
					path = original_path ,
					p_g = p_g ,
					p_t = p_t ,
					p_a = p_a ,
					p_c = p_c ,
					gc_skew = gc_sk,
					at_skew = at_sk)
				write_output(file = mito_file, 
					family = family, 
					subfamily = subfamily, 
					tribe = tribe, 
					species_name = species_name, 
					length = mito_length, 
					gc = mito_gc, 
					sex = sex_type, 
					gene_char = gene_characteristics, 
					path = original_path)
	if args.Cai:
		calculate_cai(total_gene_dict, args.GeneCode, original_path)
	os.chdir(original_path + '/' + args.MitoFile)

def read_mitogenomes(mito_file, path, file):
	print('Extracting info on %s' %(mito_file[:-4]), end = '\r')
	#print(str(i)+'/'+str(len(record['IdList'])), end='\r')
	mitogenome = SeqIO.read(mito_file,'fasta')
	os.chdir(path + '/' + args.AnnFile + '/' + file)
	for ann_file in sorted(os.listdir()):
		if ann_file.endswith('geneAnnotation.tab'):
			if ann_file[:-19] == mito_file[:-4]:
				annotation = open(ann_file,'r')
		elif ann_file.endswith('geneAnnotations.tab'):
			if ann_file[:-20] == mito_file[:-4]:
				annotation = open(ann_file,'r')
	os.chdir(path + '/' + args.MitoFile + '/' + file)
	return mitogenome, annotation

def read_taxonid(csv_file, mito_file, sex, original_path):
	samples = read_excel(os.path.join(original_path,csv_file))
	samples = np.array(samples)
	if sex == 'F':
		ncol = 4
	elif sex == 'M':
		ncol = 7
	family = ''
	subfamily = ''
	tribe = ''
	species_name = ''
	for i in range(len(samples[:,ncol])):
		if samples[i,ncol] == mito_file[:-4]:
			family = samples[i,0]
			subfamily = samples[i,1]
			tribe = samples[i,2]
			species_name = samples[i,3]
	return family, subfamily, tribe, species_name

def genome_size(mitogenome):
	return len(str(mitogenome.seq))

def mitogenome_gc_content(mitogenome):
	return GC(str(mitogenome.seq))

def reverse_complement(seq):
	seq = Seq(seq)
	return str(seq.reverse_complement())

def calculate_skew(mitogenome, bool_value = 0):
	if bool_value == 0:
		seq_in = mitogenome.seq
	elif bool_value == 1:
		seq_in = mitogenome
	n_g = str(seq_in).upper().count('G')
	n_c = str(seq_in).upper().count('C')
	n_a = str(seq_in).upper().count('A')
	n_t = str(seq_in).upper().count('T')
	p_g = n_g/len(str(seq_in)) * 100
	p_c = n_c/len(str(seq_in)) * 100
	p_a = n_a/len(str(seq_in)) * 100
	p_t = n_t/len(str(seq_in)) * 100
	n_g3 = 0
	n_c3 = 0
	n_a3 = 0
	n_t3 = 0
	for base in range(len(str(seq_in))):
		n_base = base + 1
		if n_base % 3 == 0:
			if n_base == 'G':
				n_g3 += 1
			elif n_base == 'C':
				n_c3 += 1
			elif n_base == 'A':
				n_a3 += 1
			elif n_base == 'T':
				n_t3 += 1   ## estas aqui, este ciclo nao vai funcionar mt bem pq os incrementos nao vao ter em conta todas as posiçoes 3. Ver caderno, e ver como calcular 4-fold redundant sites! Calcular ATGC e AT_skew e gc_sk disso e das posiçoes 3
		else: continue
	p_g3 = n_g3/len(str(seq_in)) * 100
	p_c3 = n_c3/len(str(seq_in)) * 100
	p_a3 = n_a3/len(str(seq_in)) * 100
	p_t3 = n_t3/len(str(seq_in)) * 100
	if (n_g + n_c) != 0 and (n_a + n_t) != 0:
		gc_sk = (n_g - n_c)/(n_g + n_c)
		at_sk = (n_a - n_t)/(n_a + n_t)
	else:
		gc_sk = 0.0
		at_sk = 0.0
	if (n_g3 + n_c3) != 0 and (n_a3 + n_t3) != 0:
		gc_sk3 = (n_g3 - n_c3)/(n_g3 + n_c3)
		at_sk3 = (n_a3 - n_t3)/(n_a3 + n_t3)
	else:
		gc_sk3 = 0.0
		at_sk3 = 0.0
	return p_g, p_c, p_a, p_t, gc_sk, at_sk, p_g3, p_c3, p_a3, p_t3, gc_sk3, at_sk3


def check_ambiguity(seq):
	'''Replacing all the IUPAC ambiguity codes into N, so GUIDANCE can work properly'''
	final_seq = ''
	for i in range(len(seq)):
		if seq[i] != 'A' and seq[i] != 'T' and seq[i] != 'G' and seq[i] != 'C' and seq[i] != 'N':
			final_seq += 'N'
		else: final_seq += seq[i]
	return final_seq


def gene_calculations(mitogenome,annotation,args,file):
	'''Gene Name: Length, GC content, Start, Stop'''
	gene_dict = {}
	gene_dictionary = {}
	current_go = ''
	pcg = ['cox1','cox3','atp6','atp8','nad4l','nad4','nad6','nad1','nad5','cob','nad2','nad3','cox2']
	mitogenome_rc = mitogenome.reverse_complement()
	mitogenome_sequece = str(mitogenome.seq)
	mitogenome_rc = str(mitogenome_rc.seq)
	mitogenome_sequece = mitogenome_sequece.upper()
	mitogenome_rc = mitogenome_rc.upper()
	for lines in annotation:			
		if '---' not in lines and ' ' not in lines:
			molecule, type_, gene, start, stop, strain = lines.split('\t')
			if gene != 'FORF' and gene != 'MORF':
				if strain == '1\n':
					p_g, p_c, p_a, p_t, gc_sk, at_sk, p_g3, p_c3, p_a3, p_t3, gc_sk3, at_sk3 = calculate_skew(mitogenome_sequece[int(start)-1:int(stop)], bool_value = 1)
					if gene in pcg:
						gene_dict[gene] = [int(stop) - int(start) + 1, GC(mitogenome_sequece[int(start)-1:int(stop)]), mitogenome_sequece[int(start)-1: int(start) + 2], mitogenome_sequece[int(stop) -3:int(stop)], p_g, p_c, p_a, p_t, gc_sk, at_sk, p_g3, p_c3, p_a3, p_t3, gc_sk3, at_sk3]
					else:
						gene_dict[gene] = [int(stop) - int(start) + 1, GC(mitogenome_sequece[int(start)-1:int(stop)]), '', '', p_g, p_c, p_a, p_t, gc_sk, at_sk,'','','','','','']
					if args.Cai:
						gene_dictionary[gene] = mitogenome_sequece[int(start)-1:int(stop)]
				elif strain == '-1\n':
					p_g, p_c, p_a, p_t, gc_sk, at_sk, p_g3, p_c3, p_a3, p_t3, gc_sk3, at_sk3 = calculate_skew(mitogenome_rc[int(start)-1:int(stop)], bool_value = 1)
					if gene in pcg:
						gene_dict[gene] = [int(stop) - int(start) + 1, GC(mitogenome_rc[int(start)-1:int(stop)]), reverse_complement(mitogenome_sequece[int(stop) - 3: int(stop)]), reverse_complement(mitogenome_sequece[int(start) -1:int(start) +2]), p_g, p_c, p_a, p_t, gc_sk, at_sk, p_g3, p_c3, p_a3, p_t3, gc_sk3, at_sk3]
					else:
						gene_dict[gene] = [int(stop) - int(start) + 1, GC(mitogenome_rc[int(start)-1:int(stop)]), '', '', p_g, p_c, p_a, p_t, gc_sk, at_sk,'','','','','','']
					if args.Cai:
						gene_dictionary[gene] = mitogenome_rc[int(start)-1:int(stop)]
				current_go = current_go + gene + ': '+ strain.strip('\n') + ', '
	return gene_dict, current_go, gene_dictionary

def gene_orders(current_go):
	HF1 = 'cox1: 1, cox3: 1, atp6: 1, trnD: 1, atp8: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnQ: -1, trnC: -1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnM: -1, nad2: -1, trnE: -1, trnS1: -1, trnS2: -1, trnA: -1, trnH: 1, nad3: 1, cox2: 1,'
	MF1 = 'cox1: 1, cox3: 1, atp6: 1, trnD: 1, atp8: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnC: -1, trnQ: -1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnE: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnM: -1, nad2: -1, trnS1: -1, trnS2: -1, trnA: -1, trnH: 1, nad3: 1, cox2: 1,'
	UF1 = 'cox1: 1, cox3: 1, atp6: 1, trnD: 1, atp8: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnC: -1, trnQ: -1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnM: -1, nad2: -1, trnE: -1, trnS1: -1, trnS2: -1, trnA: -1, trnH: 1, nad3: 1, cox2: 1,'
	UF2 = 'cox1: 1, cox3: 1, atp6: 1, trnD: 1, atp8: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnC: -1, trnQ: -1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnE: -1, trnS2: -1, trnA: -1, nad3: 1, trnM: -1, nad2: -1, trnS1: -1, trnH: 1, cox2: 1,'
	UF3 = 'cox1: 1, cox3: 1, atp6: 1, nad4l: 1, nad4: 1, nad1: -1, trnL2: -1, trnQ: -1, trnD: 1, atp8: 1, nad6: -1, trnG: -1, trnV: -1, trnI: -1, trnC: -1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnE: -1, trnS2: -1, trnA: -1, nad3: 1, trnM: -1, nad2: -1, trnS1: -1, trnH: 1, cox2: 1,'
	HM1 = 'cox1: 1, cox3: 1, atp6: 1, atp8: 1, trnD: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, trnQ: -1, nad1: -1, trnL2: -1, trnV: -1, trnH: 1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnM: -1, nad2: -1, trnE: -1, trnS1: -1, trnI: -1, trnS2: -1, trnA: -1, trnC: -1, nad3: 1, cox2: 1,'
	MM1 = 'cox1: 1, cox3: 1, atp6: 1, atp8: 1, trnD: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnC: -1, trnQ: -1, trnH: 1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnS1: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnM: -1, nad2: -1, trnE: -1, trnS2: -1, trnA: -1, nad3: 1, cox2: 1,'
	UM1 = 'cox1: 1, cox3: 1, atp6: 1, atp8: 1, trnD: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnC: -1, trnQ: -1, trnH: 1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnM: -1, nad2: -1, trnE: -1, trnS1: -1, trnS2: -1, trnA: -1, nad3: 1, cox2: 1,'
	final_go = 'Unrecognized'
	if current_go.strip() == HF1:
		final_go = 'HF1'
	elif current_go.strip() == MF1:
		final_go = 'MF1'
	elif current_go.strip() == UF1:
		final_go = 'UF1'
	elif current_go.strip() == UF2:
		final_go = 'UF2'
	elif current_go.strip() == UF3:
		final_go = 'UF3'
	elif current_go.strip() == HM1:
		final_go = 'HM1'
	elif current_go.strip() == MM1:
		final_go = 'MM1'
	elif current_go.strip() == UM1:
		final_go = 'UM1'
	return final_go

def calculate_cai(gene_dict, gen_code, path):
	reference = []
	print('\nPerforming codon analyses.. Please wait a few minutes')
	pcg = ['cox1','cox3','atp6','atp8','nad4l','nad4','nad6','nad1','nad5','cob','nad2','nad3','cox2']
	for sample in gene_dict.keys():
		print(sample)
		for gene in gene_dict[sample].keys():
			print(gene)
			cai_value = 'NA'
			if gene in pcg:
				sequence = gene_dict[sample][gene]
				sequence = check_ambiguity(sequence)
				if len(sequence) % 3 == 0:
					for other_samples in gene_dict.keys():
						if other_samples != sample and gene in gene_dict[other_samples].keys():
							if len(gene_dict[other_samples][gene]) % 3 == 0:
								reference.append(gene_dict[other_samples][gene])
					cai_value = CAI(sequence, reference = reference, genetic_code = int(gen_code))
					write_evol_characteristics(file = sample, gene = gene, cai = cai_value, path = path)

def write_evol_characteristics(file = '', gene = '', cai = 0, path = ''):
	if not os.path.exists(os.path.join(path, 'Mitogenome_Statistics/Mitogenomes_advanced.xlsx')):
		final_file = open(os.path.join(path, 'Mitogenome_Statistics/Mitogenomes_advanced.xlsx'), 'a+')
		final_file.write('Sample,Gene,CAI\n')
		final_file.close()
	final_file = open(os.path.join(path, 'Mitogenome_Statistics/Mitogenomes_advanced.xlsx'), 'a+')
	final_file.write(file + ',' + str(gene) + ',' + str(cai) + '\n')
	final_file.close()

def write_output(file = '', content = 'Gene', family = '', subfamily = '', tribe = '', species_name = '', sex = '', length = '', gc = '', gene_char = '', order = '', path = '', p_g = '', p_c = '', p_a = '', p_t = '', at_skew = '', gc_skew = ''):
	if not os.path.exists(os.path.join(path, 'Mitogenome_Statistics/Mitogenomes.xlsx')):
		final_file = open(os.path.join(path, 'Mitogenome_Statistics/Mitogenomes.xlsx'), 'a+')
		final_file.write('Sample,Content,Sex,Family,Subfamily,Tribe,Species,Gene Order,Length,%G,%C,%A,%T,GC_skew,AT_skew,GC content,Start,Stop,%G3,%C3,%A3,%T3,GC_skew3,AT_skew3\n')
		final_file.close()
	if content == 'Mitogenome':
		final_file = open(os.path.join(path, 'Mitogenome_Statistics/Mitogenomes.xlsx'), 'a+')
		final_file.write(file[:-4] + ',' + content + ',' + sex + ',' + str(family) + ',' + str(subfamily) + ',' + str(tribe) + ',' + species_name + ',' + order + ',' + str(length) + ',' + str(p_g) + ',' + str(p_c) + ',' + str(p_a) + ',' + str(p_t) + ',' + str(gc_skew) + ',' + str(at_skew) + ',' + str(gc) + ', , ' + '\n')
		final_file.close()
	elif content == 'Gene':
		final_file = open(os.path.join(path, 'Mitogenome_Statistics/Mitogenomes.xlsx'), 'a+')
		for gene in sorted(gene_char):
			final_file.write(file[:-4] + ',' + gene + ',' + sex + ',' + str(family) + ',' + str(subfamily) + ',' + str(tribe) + ',' + species_name + ',' +  order + ',' + str(gene_char[gene][0]) + ',' + str(gene_char[gene][4]) + ',' + str(gene_char[gene][5]) + ',' + str(gene_char[gene][6]) + ',' + str(gene_char[gene][7]) + ',' + str(gene_char[gene][8]) + ',' + str(gene_char[gene][9]) + ',' + str(gene_char[gene][1]) + ',' + str(gene_char[gene][2]) + ',' + str(gene_char[gene][3]) + ',' + str(gene_char[gene][10]) + ',' + str(gene_char[gene][11]) + ',' + str(gene_char[gene][12]) + ',' + str(gene_char[gene][13]) + ',' + str(gene_char[gene][14]) + ',' + str(gene_char[gene][15]) + '\n')
		final_file.close()

if __name__ == "__main__":
	args = parse()
	main(args)