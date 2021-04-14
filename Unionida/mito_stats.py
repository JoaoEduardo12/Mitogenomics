import sys, os, subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio.SeqUtils import GC
Entrez.email='joao.edu.t@hotmail.com'


def parse():
	parser = argparse.ArgumentParser(prog="gene_stats")
	parser = argparse.ArgumentParser(description='Give both the directories of where the mitogenomes are and where the annotation tables are')
	parser.add_argument('-mito', '--MitoFile', help='Give this command the directory of mitogenomes')
	parser.add_argument('-ann', '--AnnFile', help='Give this command the directory of annotations')
	args = parser.parse_args()
	return args

def main(args):
	'''Iterates over every directory and calls all statistical functions'''
	original_path = os.getcwd()
	os.mkdir('Mitogenome_Statistics')
	os.chdir(args.MitoFile)
	for file in sorted(os.listdir()):
		print(file)
		mitogenome = SeqIO.read(file,'fasta')
		os.chdir(original_path + '/' + args.AnnFile)
		for ann_name in sorted(os.listdir()):
			if ann_name.endswith('geneAnnotation.tab'):
				if ann_name[:-19] == file[:-4]:
					annotation = open(ann_name,'r')
			elif ann_name.endswith('geneAnnotations.tab'):
				if ann_name[:-20] == file[:-4]:
					annotation = open(ann_name,'r')
		mito_len = genome_size(mitogenome)
		mito_gc = mitogenome_gc_content(mitogenome)
		p_g, p_c, p_a, p_t, gc_sk, at_sk = calculate_skew(mitogenome)
		gene_characteristics, current_go = gene_calculations(mitogenome, annotation)
		mitogenome_order = gene_orders(current_go)
		write_output(file, mito_len, mito_gc, gene_characteristics, mitogenome_order, original_path)
		write_mito_excel(file, mito_len, mito_gc, mitogenome_order, p_g, p_c, p_a, p_t, gc_sk, at_sk, original_path)
		write_genelen_excel(file, gene_characteristics, original_path)
		write_genegc_excel(file, gene_characteristics, original_path)
		write_genestart_excel(file, gene_characteristics, original_path)
		write_genestop_excel(file, gene_characteristics, original_path)
		write_pg_excel(file, gene_characteristics, original_path)
		write_pc_excel(file, gene_characteristics, original_path)
		write_pa_excel(file, gene_characteristics, original_path)
		write_pt_excel(file, gene_characteristics, original_path)
		write_gcsk_excel(file, gene_characteristics, original_path)
		write_atsk_excel(file, gene_characteristics, original_path)
		os.chdir(original_path + '/' + args.MitoFile)


def genome_size(mitogenome):
	return len(str(mitogenome.seq))

def mitogenome_gc_content(mitogenome):
	return GC(str(mitogenome.seq))

def reverse_complement(seq):
	seq = Seq(seq)
	return str(seq.reverse_complement())

def calculate_skew(mitogenome, bool_value = 0):
	if bool_value == 0:
		n_g = str(mitogenome.seq).upper().count('G')
		n_c = str(mitogenome.seq).upper().count('C')
		n_a = str(mitogenome.seq).upper().count('A')
		n_t = str(mitogenome.seq).upper().count('T')
		p_g = n_g/len(str(mitogenome.seq)) * 100
		p_c = n_c/len(str(mitogenome.seq)) * 100
		p_a = n_a/len(str(mitogenome.seq)) * 100
		p_t = n_t/len(str(mitogenome.seq)) * 100
	elif bool_value == 1:
		n_g = str(mitogenome).upper().count('G')
		n_c = str(mitogenome).upper().count('C')
		n_a = str(mitogenome).upper().count('A')
		n_t = str(mitogenome).upper().count('T')
		p_g = n_g/len(str(mitogenome)) * 100
		p_c = n_c/len(str(mitogenome)) * 100
		p_a = n_a/len(str(mitogenome)) * 100
		p_t = n_t/len(str(mitogenome)) * 100
	if (n_g + n_c) != 0 and (n_a + n_t) != 0:
		gc_sk = (n_g - n_c)/(n_g + n_c)
		at_sk = (n_a - n_t)/(n_a + n_t)
	else:
		gc_sk = 0.0
		at_sk = 0.0
	return p_g, p_c, p_a, p_t, gc_sk, at_sk

def gene_calculations(mitogenome,annotation):
	'''Gene Name: Length, GC content, Start, Stop'''
	gene_dict = {}
	current_go = ''
	mitogenome_rc = mitogenome.reverse_complement()
	mitogenome_sequece = str(mitogenome.seq)
	mitogenome_rc = str(mitogenome_rc.seq)
	mitogenome_sequece = mitogenome_sequece.upper()
	mitogenome_rc = mitogenome_rc.upper()
	for lines in annotation:			
		if '---' not in lines and ' ' not in lines:
			molecule, type_, gene, start, stop, strain = lines.split('\t')
			if gene != 'FORF':
				if strain == '1\n':
					p_g, p_c, p_a, p_t, gc_sk, at_sk = calculate_skew(mitogenome_sequece[int(start):int(stop)], bool_value = 1)
					gene_dict[gene] = [int(stop) - int(start), GC(mitogenome_sequece[int(start):int(stop)]), mitogenome_sequece[int(start)-1: int(start) + 2], mitogenome_sequece[int(stop) -3:int(stop)], p_g, p_c, p_a, p_t, gc_sk, at_sk]
				elif strain == '-1\n':
					p_g, p_c, p_a, p_t, gc_sk, at_sk = calculate_skew(mitogenome_rc[int(start):int(stop)], bool_value = 1)
					gene_dict[gene] = [int(stop) - int(start), GC(mitogenome_rc[int(start):int(stop)]), reverse_complement(mitogenome_sequece[int(stop) - 3: int(stop)]), reverse_complement(mitogenome_sequece[int(start) -1:int(start) +2]), p_g, p_c, p_a, p_t, gc_sk, at_sk]
				current_go = current_go + gene + ': '+ strain.strip('\n') + ', '
	return gene_dict, current_go

def gene_orders(current_go):
	HF1 = 'cox1: 1, cox3: 1, atp6: 1, trnD: 1, atp8: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnQ: -1, trnC: -1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnM: -1, nad2: -1, trnE: -1, trnS1: -1, trnS2: -1, trnA: -1, trnH: 1, nad3: 1, cox2: 1,'
	MF1 = 'cox1: 1, cox3: 1, atp6: 1, trnD: 1, atp8: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnC: -1, trnQ: -1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnE: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnM: -1, nad2: -1, trnS1: -1, trnS2: -1, trnA: -1, trnH: 1, nad3: 1, cox2: 1,'
	UF1 = 'cox1: 1, cox3: 1, atp6: 1, trnD: 1, atp8: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnC: -1, trnQ: -1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnM: -1, nad2: -1, trnE: -1, trnS1: -1, trnS2: -1, trnA: -1, trnH: 1, nad3: 1, cox2: 1,'
	UF2 = 'cox1: 1, cox3: 1, atp6: 1, trnD: 1, atp8: 1, nad4l: 1, nad4: 1, nad6: -1, trnG: -1, nad1: -1, trnL2: -1, trnV: -1, trnI: -1, trnC: -1, trnQ: -1, nad5: 1, trnF: -1, cob: -1, trnP: -1, trnN: -1, trnL1: -1, rrnL: -1, trnY: -1, trnT: -1, trnK: -1, rrnS: -1, trnR: -1, trnW: -1, trnE: -1, trnS2: -1, trnA: -1, nad3: 1, trnM: -1, nad2: -1, trnS1: -1, trnH: 1, cox2: 1,'
	final_go = ''
	if current_go.strip() == HF1:
		final_go = 'HF1'
	if current_go.strip() == MF1:
		final_go = 'MF1'
	if current_go.strip() == UF1:
		final_go = 'UF1'
	if current_go.strip() == UF2:
		final_go = 'UF2'
	return final_go

def write_output(filename, mito_len, mito_gc, gene_characteristics, mitogenome_orders, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/mitogenome_statistics.txt')):
		final_file = open(os.path.join(original_path,'Mitogenome_Statistics/mitogenome_statistics.txt'),'a+')
		make_intro(final_file)
		final_file.write('\n\n\n================= %s ===================\n\n' % (filename))
		final_file.write('Mitogenome Length: %s' % (mito_len))
		final_file.write('\n\nMitogenome GC content: %s \n\n' % (mito_gc))
		final_file.write('\n\nGene|Length|GC content|Start & Stop Codons\n')
		for gene in gene_characteristics:
			final_file.write(gene + '|' + str(gene_characteristics[gene][0]) + '|' + str(gene_characteristics[gene][1]) + '|' + str(gene_characteristics[gene][2]) + '...' + str(gene_characteristics[gene][3]) + '\n')
		final_file.write('\n\nGene Order: %s' % (mitogenome_orders))
		final_file.close()
	else:
		final_file = open(os.path.join(original_path,'Mitogenome_Statistics/mitogenome_statistics.txt'),'a+')
		final_file.write('\n\n\n================= %s ===================\n\n' % (filename))
		final_file.write('Mitogenome Length: %s' % (mito_len))
		final_file.write('\n\nMitogenome GC content: %s \n\n' % (mito_gc))
		final_file.write('\n\nGene|Length|GC content|Start & Stop Codons\n')
		for gene in gene_characteristics:
			final_file.write(gene + '|' + str(gene_characteristics[gene][0]) + '|' + str(gene_characteristics[gene][1]) + '|' + str(gene_characteristics[gene][2]) + '...' + str(gene_characteristics[gene][3]) + '\n')
		final_file.write('\n\nGene Order: %s' % (mitogenome_orders))

def make_intro(file):
	first = "=" * 50
	second = 'Mitogenome  Analyses'
	third = "=" * 15
	file.write(first + '\n' + first + '\n' + third + second + third + '\n' + first + '\n' + first + '\n' + '\n')

def write_mito_excel(filename, mito_len, mito_gc, mitogenome_order,p_g, p_c, p_a, p_t, gc_sk, at_sk, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/mitogenome_statistics.xlsx')):
		mito_csv = open(os.path.join(original_path,'Mitogenome_Statistics/mitogenome_statistics.xlsx'),'a+')
		mito_csv.write('Sequence Code,Gene Order,Mitogenome Length,%G,%C,%A,%T,GC_skew,AT_skew,GC content\n')
		mito_csv.write(filename + ',' + mitogenome_order + ',' + str(mito_len) + ',' + str(p_g) + ',' + str(p_c) + ',' + str(p_a) + ',' + str(p_t) + ',' + str(gc_sk) + ',' + str(at_sk) + ',' + str(mito_gc) + '\n')
		mito_csv.close()
	else:
		mito_csv = open(os.path.join(original_path,'Mitogenome_Statistics/mitogenome_statistics.xlsx'),'a+')
		mito_csv.write(filename + ',' + mitogenome_order + ',' + str(mito_len) + ',' + str(p_g) + ',' + str(p_c) + ',' + str(p_a) + ',' + str(p_t) + ',' + str(gc_sk) + ',' + str(at_sk) + ',' + str(mito_gc) + '\n')
		mito_csv.close()

def write_genelen_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_length.xlsx')):
		gene_len_csv = open(os.path.join(original_path, 'Mitogenome_Statistics/gene_length.xlsx'),'a+')
		gene_len_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_len_csv.write(gene + ',')
		gene_len_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_len_csv.write(',' + str(gene_characteristics[gene][0]))
		gene_len_csv.close()
	else:
		gene_len_csv = open(os.path.join(original_path, 'Mitogenome_Statistics/gene_length.xlsx'),'a+')
		gene_len_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_len_csv.write(',' + str(gene_characteristics[gene][0]))
		gene_len_csv.close()

def write_genegc_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_gc.xlsx')):
		gene_gc_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_gc.xlsx'),'a+')
		gene_gc_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_gc_csv.write(gene + ',')
		gene_gc_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_gc_csv.write(',' + str(gene_characteristics[gene][1]))
		gene_gc_csv.close()
	else:
		gene_gc_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_gc.xlsx'),'a+')
		gene_gc_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_gc_csv.write(',' + str(gene_characteristics[gene][1]))
		gene_gc_csv.close()

def write_genestart_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_start.xlsx')):
		gene_start_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_start.xlsx'),'a+')
		gene_start_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_start_csv.write(gene + ',')
		gene_start_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_start_csv.write(',' + str(gene_characteristics[gene][2]))
		gene_start_csv.close()
	else:
		gene_start_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_start.xlsx'),'a+')
		gene_start_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_start_csv.write(',' + str(gene_characteristics[gene][2]))
		gene_start_csv.close()

def write_genestop_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_stop.xlsx')):
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_stop.xlsx'),'a+')
		gene_stop_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(gene + ',')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][3]))
		gene_stop_csv.close()
	else:
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_stop.xlsx'),'a+')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][3]))
		gene_stop_csv.close()

def write_pg_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_G.xlsx')):
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_G.xlsx'),'a+')
		gene_stop_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(gene + ',')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][4]))
		gene_stop_csv.close()
	else:
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_G.xlsx'),'a+')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][4]))
		gene_stop_csv.close()

def write_pc_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_C.xlsx')):
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_C.xlsx'),'a+')
		gene_stop_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(gene + ',')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][5]))
		gene_stop_csv.close()
	else:
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_C.xlsx'),'a+')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][5]))
		gene_stop_csv.close()

def write_pa_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_A.xlsx')):
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_A.xlsx'),'a+')
		gene_stop_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(gene + ',')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][6]))
		gene_stop_csv.close()
	else:
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_A.xlsx'),'a+')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][6]))
		gene_stop_csv.close()

def write_pt_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_T.xlsx')):
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_T.xlsx'),'a+')
		gene_stop_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(gene + ',')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][7]))
		gene_stop_csv.close()
	else:
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_T.xlsx'),'a+')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][7]))
		gene_stop_csv.close()

def write_gcsk_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_GCskew.xlsx')):
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_GCskew.xlsx'),'a+')
		gene_stop_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(gene + ',')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][8]))
		gene_stop_csv.close()
	else:
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_GCskew.xlsx'),'a+')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][8]))
		gene_stop_csv.close()

def write_atsk_excel(filename, gene_characteristics, original_path):
	if not os.path.exists(os.path.join(original_path, 'Mitogenome_Statistics/gene_ATskew.xlsx')):
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_ATskew.xlsx'),'a+')
		gene_stop_csv.write('Sequence Code/Genes,')
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(gene + ',')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][9]))
		gene_stop_csv.close()
	else:
		gene_stop_csv = open(os.path.join(original_path,'Mitogenome_Statistics/gene_ATskew.xlsx'),'a+')
		gene_stop_csv.write('\n' + filename)
		for gene in sorted(gene_characteristics):
			gene_stop_csv.write(',' + str(gene_characteristics[gene][9]))
		gene_stop_csv.close()

if __name__ == "__main__":
	args = parse()
	main(args)