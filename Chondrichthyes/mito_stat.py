import sys, os, subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio import Entrez
from Bio.GenBank import Record
from CAI import RSCU
from CAI import CAI
from collections import OrderedDict
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Phylo.PAML import codeml
Entrez.email='joao.edu.t@hotmail.com'
api_key = '105bc09b7b34973aa305f31be73aecbc4109'

def parse_args():
	"""Parsing arguments via terminal"""
	parser = argparse.ArgumentParser(prog="mitogenome_statistics")
	parser = argparse.ArgumentParser(description='Given a search term in NCBI, reads all available sequences and extracts multiple statistics')
	parser.add_argument('-search', '--Search', help='Give this command a file with the search term in NCBI (right side of the page in "Search Details, when you query an entry")')
	parser.add_argument('-taxonomy', '--Tax', help = 'Give this command the highest recent common taxonomy classification within search (optional)')
	parser.add_argument('-genetic_code', '--GeneCode', help = 'Give this command the genetic code of the input sequences. Defaults to 11')
	parser.add_argument('-download', "--DownloadDir", help = '(Optional) Give this command a directory name where the complete mitogenomes will be saved as well as each gene fasta')
	parser.add_argument('-clean', "--Clean", help = '(Optional) Give this command the argument to remove sequences which have at least one gene with gaps')
	parser.add_argument('-remove_accessions', "--Accessions",
						help='Give this command a file with accession numbers marked with hashtag beforehand (anywhere)')
	args = parser.parse_args()
	return args

def main(args):
	"""Uses the NCBI search utilities to retrieve sequences and genes inside"""
	if args.DownloadDir:
		os.mkdir(args.DownloadDir)
	unverified, no_features, no_genes = [], [], []
	print('\nSearching in NCBI... Hold tight!\n')
	handle = Entrez.esearch(db='nucleotide', retmax=100000, term= args.Search, idtype="acc")
	record = Entrez.read(handle)
	print('\nThere are %s results found\n' % str(len(record['IdList'])))
	handle.close()
	print('\nExtracting features.. this may take a while\n')
	i = 0
	full_gene_dict = {}
	for rec in record['IdList']:
		i += 1
		print(' ' + str(i)+'/'+str(len(record['IdList'])), end='\r')
		unverified, no_features, no_genes, full_gene_dict = extract_features(rec, unverified, no_features, no_genes, full_gene_dict, args)
	calculate_cai(full_gene_dict, args.GeneCode)
	calculate_dnds(full_gene_dict)
	write_log(unverified = unverified, no_features = no_features, missing_gene = no_genes)

def extract_features(record, unverified, no_features, no_genes, full_gene_dict, args):
	clean = False
	if args.Accessions:
		records = filter_accessions(args.Accessions)
	else: records = []
	mitogenome = Entrez.efetch(db = "nucleotide", id = record, rettype = 'fasta', retmode = 'text')
	mitogenome_sequence = mitogenome.read()
	mitogenome_characteristics = mitogenome_calculations(mitogenome_sequence)
	#if args.DownloadDir:
	#	file = open(os.path.join(args.DownloadDir, gene) + ".fn", 'a+')
	#	file.write(">%s\n%s\n" % (name, mitogenome_sequence))
	#	file.close()
	gene_dict = {}
	state = ''
	with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=record) as handle:
		for rec in SeqIO.parse(handle,'gb'):
			if rec.id not in records:
				taxonomy_names = ';'.join(rec.annotations['taxonomy'])
				if 'UNVERIFIED' in rec.description:
					unverified.append(rec.description)
					state = 'unverified'
				elif 'partial' in rec.description:
					state = 'partial'
				else:
					state = 'complete'
				if rec.features:
					for feature in rec.features:
						if feature.type == 'CDS':
							if 'gene' in feature.qualifiers.keys():
								gene_sequence = feature.location.extract(rec).seq
								if args.Clean:
									clean = filter_sequences(gene_sequence, rec.annotations['organism'], full_gene_dict, feature.qualifiers['gene'][0].strip(''))
									if clean:
										gene_characteristics, full_gene, gene = gene_calculations(gene_dictionary = gene_dict,
											gene_name = feature.qualifiers['gene'][0].strip(''),
											sequence = gene_sequence, pcg = 1, args = args, accession = rec.id)
										full_gene_dict[(rec.annotations['organism'],rec.id,taxonomy_names.split(args.Tax + ';')[1],state, gene)] = full_gene
								else:
									gene_characteristics, full_gene, gene = gene_calculations(gene_dictionary=gene_dict,
																							  gene_name=
																							  feature.qualifiers['gene'][
																								  0].strip(''),
																							  sequence=gene_sequence, pcg=1,
																							  args=args, accession=rec.id)
									full_gene_dict[(
									rec.annotations['organism'], rec.id, taxonomy_names.split(args.Tax + ';')[1], state,
									gene)] = full_gene
							else:
								gene_sequence = feature.location.extract(rec).seq
								if args.Clean:
									clean = filter_sequences(gene_sequence, rec.annotations['organism'], full_gene_dict, feature.qualifiers['product'][0].strip(''))
									if clean:
										gene_characteristics, full_gene, gene = gene_calculations(gene_dictionary = gene_dict,
											gene_name = feature.qualifiers['product'][0].strip(''),
											sequence = gene_sequence, pcg = 1, args = args, accession = rec.id)
										full_gene_dict[(rec.annotations['organism'],rec.id,taxonomy_names.split(args.Tax + ';')[1],state, gene)] = full_gene
									else:
										gene_characteristics, full_gene, gene = gene_calculations(gene_dictionary=gene_dict,
																								  gene_name=
																								  feature.qualifiers[
																									  'product'][0].strip(
																									  ''),
																								  sequence=gene_sequence,
																								  pcg=1, args=args,
																								  accession=rec.id)
										full_gene_dict[(
										rec.annotations['organism'], rec.id, taxonomy_names.split(args.Tax + ';')[1], state,
										gene)] = full_gene

						elif feature.type == 'tRNA':
							gene_sequence = feature.location.extract(rec).seq
							if args.Clean:
								if 'product' in feature.qualifiers.keys():
									clean = filter_sequences(gene_sequence, rec.annotations['organism'], full_gene_dict, feature.qualifiers['product'][0].strip(''))
									if clean:
										gene_characteristics, full_gene, gene = gene_calculations(gene_dictionary = gene_dict,
											gene_name = feature.qualifiers['product'][0].strip(''),
											sequence = gene_sequence, args = args, accession = rec.id)
							else:
								if 'product' in feature.qualifiers.keys():
									gene_characteristics, full_gene, gene = gene_calculations(gene_dictionary=gene_dict,
																						  gene_name=
																						  feature.qualifiers['product'][
																							  0].strip(''),
																						  sequence=gene_sequence, args=args,
																						  accession=rec.id)

								#full_gene_dict[(rec.annotations['organism'],rec.id,taxonomy_names.split(args.Tax + ';')[1],state, gene)] = full_gene
							#else:
							#	gene_characteristics = gene_calculations(gene_dictionary = gene_dict,
							#		gene_name = feature.qualifiers['gene'][0].strip(''),
							#		sequence = gene_sequence)
						elif feature.type == 'rRNA':
							gene_sequence = feature.location.extract(rec).seq
							if args.Clean:
								clean = filter_sequences(gene_sequence, rec.annotations['organism'], full_gene_dict, feature.qualifiers['product'][0].strip(''))
								if clean:
									if 'product' in feature.qualifiers.keys():
										gene_characteristics, full_gene, gene = gene_calculations(gene_dictionary = gene_dict,
											gene_name = feature.qualifiers['product'][0].strip(''),
											sequence = gene_sequence, args = args, accession = rec.id)
							else:
								if 'product' in feature.qualifiers.keys():
									gene_characteristics, full_gene, gene = gene_calculations(gene_dictionary=gene_dict,
																							  gene_name=
																							  feature.qualifiers['product'][
																								  0].strip(''),
																							  sequence=gene_sequence,
																							  args=args, accession=rec.id)
								#full_gene_dict[(rec.annotations['organism'],rec.id,taxonomy_names.split(args.Tax + ';')[1],state, gene)] = full_gene
							#else:
							#	gene_characteristics = gene_calculations(gene_dictionary = gene_dict,
							#		gene_name = feature.qualifiers['gene'][0].strip(''),
							#		sequence = gene_sequence)
						elif feature.type == 'source':
							gene_characteristics = []
				else:
					no_features.append(rec.description)
					gene_characteristics = []
				assembly_method = ''
				sequencing_technology = ''
				if 'structured_comment' in rec.annotations:
					if 'Assembly-Data' in rec.annotations['structured_comment'].keys():
						if 'Assembly Method' in rec.annotations['structured_comment']['Assembly-Data'].keys():
							assembly_method = rec.annotations['structured_comment']['Assembly-Data']['Assembly Method']
						if 'Sequencing Technology' in rec.annotations['structured_comment']['Assembly-Data'].keys():
							sequencing_technology = rec.annotations['structured_comment']['Assembly-Data']['Sequencing Technology']
				if args.Tax:
					taxonomy_names = taxonomy_names.split(args.Tax + ';')[1]
				if gene_characteristics != []:
					gene_check(gene_characteristics,no_genes,rec.annotations['organism'],taxonomy_names)
				if args.Clean:
					if clean:
						write_output(name=rec.annotations['organism'],
									 acc_number=rec.id,
									 family=taxonomy_names,
									 mito_char=mitogenome_characteristics,
									 gene_char=gene_characteristics,
									 state=state,
									 assembly=assembly_method,
									 sequencing=sequencing_technology)
				else: write_output(name = rec.annotations['organism'],
					acc_number = rec.id,
					family = taxonomy_names,
					mito_char = mitogenome_characteristics,
					gene_char = gene_characteristics,
					state = state,
					assembly = assembly_method,
					sequencing = sequencing_technology)
	return unverified, no_features, no_genes, full_gene_dict

def filter_accessions(accessions):
	records = []
	with open(accessions) as f:
		lines = f.readlines()
		f.close()
	for line in lines:
		if '#' in line:
			records.append(line.strip('\n')[1:])
	return records

def filter_sequences(seq, name, group, gene):
	clean = True
	#print('\nThis is the sequence: ', seq)
	#print('\nThis is the name: ', name)
	#print('\nThis is the group: ',group)
	#print('\nThis is the gene: ', gene)
	for key in group:
		if (key[0] == name and key[4] == gene) or seq.count('N') > 7:
			print('\n{} is no good, removing it from the dataset' .format(key[1]))
			clean = False
	#if (name in group.keys(0) and gene in group.keys(0)) or seq.count('N') > 7:
	#	print('\n\n\n\n\n\nWE DID IT BOYS\n\n\n\n\n\n\n\n')
	#	clean = False
	return clean

def mitogenome_calculations(sequence):
	'''1- Genome Size; 2- Genome GC content; 3-G; 4-C; 5-A; 6-T; 7- GC skew; 8- AT skew'''
	mitogenome_characteristics = []
	mitogenome_characteristics.append(len(sequence))
	mitogenome_characteristics.append(GC(sequence))
	n_g = str(sequence).upper().count('G')
	n_c = str(sequence).upper().count('C')
	n_a = str(sequence).upper().count('A')
	n_t = str(sequence).upper().count('T')
	p_g = n_g/len(str(sequence)) * 100
	mitogenome_characteristics.append(p_g)
	p_c = n_c/len(str(sequence)) * 100
	mitogenome_characteristics.append(p_c)
	p_a = n_a/len(str(sequence)) * 100
	mitogenome_characteristics.append(p_a)
	p_t = n_t/len(str(sequence)) * 100
	mitogenome_characteristics.append(p_t)
	if (n_g + n_c) != 0 and (n_a + n_t) != 0:
		gc_sk = (n_g - n_c)/(n_g + n_c)
		at_sk = (n_a - n_t)/(n_a + n_t)
	else: gc_sk, at_sk = 0, 0
	mitogenome_characteristics.append(gc_sk)
	mitogenome_characteristics.append(at_sk)
	return mitogenome_characteristics

def gene_calculations(gene_dictionary = {}, gene_name = '', sequence = '', pcg = 0, args = "", accession = ""):
	gene_name = gene_correction(accession, sequence, gene_name, args)
	full_gene = ''
	n_g = str(sequence).upper().count('G')
	n_c = str(sequence).upper().count('C')
	n_a = str(sequence).upper().count('A')
	n_t = str(sequence).upper().count('T')
	p_g = n_g/len(str(sequence)) * 100
	p_c = n_c/len(str(sequence)) * 100
	p_a = n_a/len(str(sequence)) * 100
	p_t = n_t/len(str(sequence)) * 100
	if (n_g + n_c) != 0 and (n_a + n_t) != 0:
		gc_sk = (n_g - n_c)/(n_g + n_c)
		at_sk = (n_a - n_t)/(n_a + n_t)
	else: gc_sk, at_sk = 0, 0
	if pcg == 1:
		codon_list = [str(sequence)[base:base+3] for base in range(0, len(str(sequence)), 3)]
		p_g1 = ''.join([base[0] for base in codon_list if len(base) == 3]).upper().count('G')/len(codon_list) * 100
		p_c1 = ''.join([base[0] for base in codon_list if len(base) == 3]).upper().count('C')/len(codon_list) * 100
		p_a1 = ''.join([base[0] for base in codon_list if len(base) == 3]).upper().count('A')/len(codon_list) * 100
		p_t1 = ''.join([base[0] for base in codon_list if len(base) == 3]).upper().count('T')/len(codon_list) * 100
		p_g2 = ''.join([base[1] for base in codon_list if len(base) == 3]).upper().count('G')/len(codon_list) * 100
		p_c2 = ''.join([base[1] for base in codon_list if len(base) == 3]).upper().count('C')/len(codon_list) * 100
		p_a2 = ''.join([base[1] for base in codon_list if len(base) == 3]).upper().count('A')/len(codon_list) * 100
		p_t2 = ''.join([base[1] for base in codon_list if len(base) == 3]).upper().count('T')/len(codon_list) * 100
		p_g3 = ''.join([base[2] for base in codon_list if len(base) == 3]).upper().count('G')/len(codon_list) * 100
		p_c3 = ''.join([base[2] for base in codon_list if len(base) == 3]).upper().count('C')/len(codon_list) * 100
		p_a3 = ''.join([base[2] for base in codon_list if len(base) == 3]).upper().count('A')/len(codon_list) * 100
		p_t3 = ''.join([base[2] for base in codon_list if len(base) == 3]).upper().count('T')/len(codon_list) * 100
		gene_dictionary[gene_name] = [len(sequence), GC(sequence), str(p_g), str(p_c), str(p_a), str(p_t), str(gc_sk), str(at_sk), sequence[0:3], sequence[-3::], str(p_g1), str(p_c1), str(p_a1), str(p_t1), str(p_g2), str(p_c2), str(p_a2), str(p_t2), str(p_g3), str(p_c3), str(p_a3), str(p_t3)]
		full_gene = str(sequence).upper()
	else:
		if gene_name in gene_dictionary.keys():
			gene_name += '2'
		gene_dictionary[gene_name] = [len(sequence), GC(sequence), str(p_g), str(p_c), str(p_a), str(p_t), str(gc_sk), str(at_sk), '', '', '', '', '', '', '', '', '', '', '', '', '', '']
	return gene_dictionary, full_gene, gene_name

def gene_correction(name, sequence, gene, args):
	ATP6_list = ['atp6','ATP6','Atp6','Atp 6','ATP 6','ATPase 6','MTATP6','MTATP 6','atp6','ATPase 6','ATPase6','ATPase subunit 6','ATP synthase F0 subunit 6','ATP-6','atp-6','ATPase-6','atpase6','atpase 6','atpase-6','atp vi','ATP-vi','atp-vi','atpvi','atpasevi','atpase vi','Atpase VI','Atpase vi','ATP VI']
	ATP8_list = ['atp8','Atp8','Atp 8','ATPase 8','ATP 8','MTATP8','atp8','MTATP 8','ATP8','ATPase 8','ATPase8','ATP synthase F0 subunit 8','ATPase subunit 8','ATP-8','atp-8','ATPase-8','atpase8','atpase 8','atpase-8','atp viii','ATP-viii','atp-viii','atpviii','atpaseviii','atpase viii','Atpase VIII','Atpase viii','ATP VIII']
	COI_list = ['cox1','COI','CO I','cox1','COX1','COXI','COX I','COX-I','CO1','cytochrome oxidase subunit I','cytochrome c oxidase subunit I','(cytochrome oxidase subunit I)','COX-1','CO-I','CO-1','CoI','coI','co-I','co-1','co 1','Co 1','Co I','Co1','CoI','coi','co-i','co i','cox 1','cox-1','cox i','cox-i']
	COII_list = ['cox2','COII','CO II','cox2','COX2','COXII','COX II','COX-II','CO2','(cytochrome oxidase subunit II)','cytochrome c oxidase subunit II','cytochrome oxidase subunit II','COX-2','CO-II','CO-2','CoII','co-II','co-2','co 2','Co 2','Co II','Co2','CoII','coii','co-ii','co ii','cox 2','cox-2','cox ii','cox-ii']
	COIII_list = ['cox3','COIII','CO III','cox3','COX3','COX III','COX-III','COXIII','CO3','(cytochrome oxidase subunit III)','cytochrome c oxidase subunit III','cytochrome oxidase subunit III','COX-3','CO-III','CO-3','CoIII','co-III','co-3','co 3','Co 3','Co III','Co3','CoIII','coiii','co-iii','co iii','cox 3','cox-3','cox iii','cox-iii']
	CYTB_list = ['cob','CYB','CY B','cyb','cy b','CYTB','CYT B','cob','cyt b','Cyt b','cytochrome_b','cytochrome b','cytb','Cytb','CYTb','CTYB','cyt-b','CYT-B','CTY-B','cytB','cyt b','CYT b','CB','Cb','C b','C B','cb','c b','c B','cB']
	NAD1_list = ['ND1','nad1','NAD1','nadh1','NAD 1','NADH 1','NADH-1','NADH1','NADH dehydrogenase subunit 1','NAD-1','nadh-1','nad-1','nadh 1','nad 1','Nad1','Nadh1','Nad-1','Nadh-1','Nad 1','Nadh 1','NADH-I','NADH I','NAD-I','NAD I','Nad-I','Nad I','nad-I','nad I','nad-i','nadi','nad i','Nad i','Nad-i','Nadh i','Nadhi','Nadh-i','NADH i','NADH-i','NADHi','NADHI','NADI','NADi','NAD i','Nad i','nad i']
	NAD2_list = ['ND2','nad2','NAD2','NAD 2','nadh2','NADH-2','NADH 2','NADH2','NADH dehydrogenase subunit 2','NAD-2','nadh-2','nad-2','nadh 2','nad 2','Nad2','Nadh2','Nad-2','Nadh-2','Nad 2','Nadh 2','NADH-II','NADH II','NAD-II','NAD II','Nad-II','Nad II','nad-II','nad II','nad-ii','nadii','nad ii','Nad ii','Nad-ii','Nadh ii','Nadhii','Nadh-ii','NADH ii','NADH-ii','NADHii','NADHII','NADII','NADii','NAD ii','Nad ii','nad ii']
	NAD3_list = ['ND3','nad3','NAD3','NAD 3','NADH-3','nadh3','NADH 3','NADH3','NADH dehydrogenase subunit 3','NAD-3','nadh-3','nad-3','nadh 3','nad 3','Nad3','Nadh3','Nad-3','Nadh-3','Nad 3','Nadh 3','NADH-III','NADH III','NAD-III','NAD III','Nad-III','Nad III','nad-III','nad III','nad-iii','nadiii','nad iii','Nad iii','Nad-iii','Nadh iii','Nadhiii','Nadh-iii','NADH iii','NADH-iii','NADHiii','NADHIII','NADIII','NADiii','NAD iii','Nad iii','nad iii']
	NAD4_list = ['ND4','nad4','NAD4','NAD 4','nadh4','NADH-4','NADH 4','NADH4','NADH dehydrogenase subunit 4','NAD-4','nadh-4','nad-4','nadh 4','nad 4','Nad4','Nadh4','Nad-4','Nadh-4','Nad 4','Nadh 4','NADH-IV','NADH IV','NAD-IV','NAD IV','Nad-IV','Nad IV','nad-IV','nad IV','nad-iv','nadiv','nad iv','Nad iv','Nad-iv','Nadh iv','Nadhiv','Nadh-iv','NADH iv','NADH-iv','NADHiv','NADHIV','NADIV','NADiv','NAD iv','Nad iv','nad iv']
	NAD4L_list = ['ND4L','nad4l','NAD4L','NAD 4L','nadh4l','nadh4L','NADH 4L','NADH-4L','NADH4L','NADH dehydrogenase subunit 4L','NAD-4L','nadh-4L','nad-4L','nadh 4L','nad 4L','Nad4L','Nadh4L','Nad-4L','Nadh-4L','Nad 4L','Nadh 4L','NADH-IVL','NADH IVL','NAD-IVL','NAD IVL','Nad-IVL','Nad IVL','nad-IVL','nad IVL','nad-ivl','nadivl','nad ivl','Nad ivl','Nad-ivl','Nadh ivl','Nadhivl','Nadh-ivl','NADH ivl','NADH-ivl','NADHivl','NADHIVL','NADIVL','NADivl','NAD ivl','Nad ivl','nad ivl']
	NAD5_list = ['ND5','nad5','NAD5','NAD 5','NADH-5','nadh5','NADH 5','NADH5','NADH dehydrogenase subunit 5','NAD-5','nadh-5','nad-5','nadh 5','nad 5','Nad5','Nadh5','Nad-5','Nadh-5','Nad 5','Nadh 5','NADH-V','NADH V','NAD-V','NAD V','Nad-V','Nad V','nad-v','nad v','nad-v','nadv','nad v','Nad v','Nad-v','Nadh v','Nadhv','Nadh-v','NADH v','NADH-v','NADHv','NADHV','NADV','NADv','NAD v','Nad v','nad v']
	NAD6_list = ['ND6','nad6','NAD6','NAD 6','NADH-6','nadh6','NADH 6','NADH6','NADH dehydrogenase subunit 6','NAD-6','nadh-6','nad-6','nadh 6','nad 6','Nad6','Nadh6','Nad-6','Nadh-6','Nad 6','Nadh 6','NADH-VI','NADH VI','NAD-VI','NAD VI','Nad-VI','Nad VI','nad-VI','nad VI','nad-vi','nadvi','nad vi','Nad vi','Nad-vi','Nadh vi','Nadhvi','Nadh-vi','NADH vi','NADH-vi','NADHvi','NADHVI','NADVI','NADvi','NAD vi','Nad vi','nad vi']
	orfs = ['f-orf','F-ORF','F-orf','F-Orf','FORF','forf','Forf','HORF','h-orf','H-Orf','Horf','horf','H-orf','orf','Orf','ORF']
	rnaL = ['rrnL','16S','16S ribosomal RNA','16S Ribosomal RNA','l6S ribosomal RNA','l-rRNA','L-rRNA','l rRNA','L rRNA','l-rrna','l rrna','L RRNA','lrrna','l rRNA','lrRNA']
	rnaS = ['rrnS','12S','12S ribosomal RNA','12S Ribosomal RNA','l2S ribosomal RNA','s-rRNA','S-rRNA','s rRNA','S rRNA','s-rrna','s rrna','S RRNA','srrna','s rRNA','srRNA']
	d_loop = ['D-loop','D-loop','control region','d-loop','D loop','D Loop','D-Loop','D LOOP','D-LOOP','CONTROL REGION']
	tRNAs = ['tRNA-Phe','tRNA-Val','tRNA-Leu','tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Ala','tRNA-Asn','tRNA-Cys','tRNA-Tyr','tRNA-Ser','tRNA-Asp','tRNA-Lys','tRNA-Gly','tRNA-Arg','tRNA-His','tRNA-Ser','tRNA-Leu','tRNA-Glu','tRNA-Thr','tRNA-Pro']
	if gene in ATP6_list: #no case-switch in python...
		gene = ATP6_list[0]
	elif gene in ATP8_list:
		gene = ATP8_list[0]
	elif gene in COI_list:
		gene = COI_list[0]
	elif gene in COII_list:
		gene = COII_list[0]
	elif gene in COIII_list:
		gene = COIII_list[0]
	elif gene in CYTB_list:
		gene = CYTB_list[0]
	elif gene in NAD1_list:
		gene = NAD1_list[1]
	elif gene in NAD2_list:
		gene = NAD2_list[1]
	elif gene in NAD3_list:
		gene = NAD3_list[1]
	elif gene in NAD4_list:
		gene = NAD4_list[1]
	elif gene in NAD4L_list:
		gene = NAD4L_list[1]
	elif gene in NAD5_list:
		gene = NAD5_list[1]
	elif gene in NAD6_list:
		gene = NAD6_list[1]
	elif gene in rnaL:
		gene = rnaL[0]
	elif gene in rnaS:
		gene = rnaS[0]
	elif gene in orfs:
		gene = orfs[0]
	elif gene in d_loop:
		gene = d_loop[0]
	elif gene in tRNAs:
		gene = gene
	else:
		sys.exit("There seems to be a gene read (%s) that does not match the following: ATP8, ATP6, COI, COII, COIII, CYTB, ND1, ND2, ND3, ND4, ND4L, ND5, ND6" % gene)
	if args.DownloadDir:
		file = open(os.path.join(args.DownloadDir, gene) + ".fn", 'a+')
		file.write(">%s\n%s\n" % (name, sequence))
		file.close()
	return gene

def gene_check(gene_characteristics, no_genes, sample, tax):
	complete_list = ['atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6','rrnL','rrnS','tRNA-Phe','tRNA-Val','tRNA-Leu','tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Ala','tRNA-Asn','tRNA-Cys','tRNA-Tyr','tRNA-Ser','tRNA-Asp','tRNA-Lys','tRNA-Gly','tRNA-Arg','tRNA-His','tRNA-Ser','tRNA-Leu','tRNA-Glu','tRNA-Thr','tRNA-Pro','tRNA-Ser2','tRNA-Leu2']
	for gene in gene_characteristics.keys():
		if not gene in complete_list:
			gene_report = 'The gene {} is not present in {} (Taxonomy: {})' .format(gene,sample,tax)
			no_genes.append(gene_report)
		else: continue

def calculate_cai(gene_dict, gen_code):
	reference = []
	print('\nPerforming codon analyses.. Please wait a few minutes!')
	pcg = ['cox1','cox3','atp6','atp8','nad4l','nad4','nad6','nad1','nad5','cob','nad2','nad3','cox2']
	gene_dict = filter_ambiguous(gene_dict,pcg)
	for tupl in gene_dict.keys():
		for gene in pcg:
			if gene in tupl:
				cai_value = 'NA'
				sequence = gene_dict[tupl]
				if len(sequence) % 3 == 0:
					rscu_list = []
					for other_tpl in gene_dict.keys():
						if other_tpl != tupl:
							if len(gene_dict[other_tpl]) % 3 == 0:
								reference.append(gene_dict[other_tpl])
					cai_value = CAI(sequence, reference = reference, genetic_code = int(gen_code))
					rscu_list.append(sequence)
					rscu_values = RSCU(rscu_list, genetic_code = int(gen_code))
					write_evol_characteristics(file = tupl, cai = cai_value, rscu = rscu_values)

def filter_ambiguous(gene_dict, pcg):
	nucleotides = ['A','T','G','C']
	for tupl in gene_dict.keys():
		for gene in pcg:
			if gene in tupl:
				sequence = gene_dict[tupl]
				filt_sequence = ''
				for codon in [sequence[i:i+3] for i in range(0, len(sequence), 3)]:
					if len(codon) == 3:
						if codon[0] in nucleotides and codon[1] in nucleotides and codon[2] in nucleotides:
							filt_sequence += codon
						else: continue
				gene_dict[tupl] = filt_sequence
	return gene_dict

def calculate_dnds(all_gene_dict):
	os.mkdir('protein_coding')
	os.mkdir('amino_acids')
	os.mkdir('aa_alignments')
	os.mkdir('codon_alignments')
	os.mkdir('KaKs')
	curr_dir = os.getcwd()
	pcg = ['cox1','cox3','atp6','atp8','nad4l','nad4','nad6','nad1','nad5','cob','nad2','nad3','cox2']
	for tupl in all_gene_dict.keys():
		for gene in pcg:
			if gene in tupl:
				gene_file = open('protein_coding/{}.fna' .format(gene), 'a+')
				gene_file.write('>{}\n' .format(tupl[1]))
				gene_file.write(all_gene_dict[tupl] + '\n')
				gene_file.close()
				aa_file = open('amino_acids/{}.faa' .format(gene), 'a+')
				aa_file.write('>{}\n' .format(tupl[1]))
				aa_file.write(str(Seq(all_gene_dict[tupl]).translate(table = args.GeneCode)) + '\n')
				aa_file.close()
	for gene in pcg:
		os.mkdir('KaKs/{}' .format(gene))
		muscle_cline = MuscleCommandline(input = 'amino_acids/{}.faa' .format(gene)) #out = 'aa_alignments/{}.aln.faa' .format(gene)
		stdout, stderr = muscle_cline()
		with open('aa_alignments/{}.aln.faa' .format(gene), 'a+') as handle:
			handle.write(stdout)
		align = AlignIO.read('aa_alignments/{}.aln.faa' .format(gene), 'fasta')
		ctl_file = open('codeml.ctl', 'a+')
		ctl_file.write('seqfile = {}.pal2nal\noutfile = {}_codeml.txt\nnoisy = 0\nverbose = 1\nrunmode = -2\nseqtype = 1\nCodonFreq = 2\nmodel = 1\nNSsites = 1\nicode = {}\nfix_kappa = 1\nkappa = 1\nfix_omega = 0\nomega = 0.5\ncleandata = 1' .format(gene,gene,str(int(args.GeneCode) -1)))
		ctl_file.close()
		#os.system('/home/edu/Desktop/Softwares/guidance.v2.02/www/Guidance/guidance.pl --seqFile protein.fas --msaProgram MAFFT --seqType aa --outDir /somewhere/protein.guidance')
		os.system('dnds/pal2nal.pl aa_alignments/{}.aln.faa protein_coding/{}.fna -output paml -codontable {} -nogap > {}.pal2nal' .format(gene,gene,args.GeneCode,gene))
		os.system('codeml')
		os.system('python dnds/parse_codeml_output.py {}_codeml.txt' .format(gene))
		make_dnds_csv('{}_codeml.txt' .format(gene), all_gene_dict, gene)
		os.system('rm codeml.ctl')
		os.system('mv {}_dnds.csv KaKs/{}' .format(gene,gene))
		os.system('mv *.pal2nal codon_alignments')
		os.system('mv *_codeml.txt KaKs/{}' .format(gene))
		os.system('mv rst KaKs/{}' .format(gene))
		os.system('mv rst1 KaKs/{}' .format(gene))
		os.system('mv rub KaKs/{}' .format(gene))
		os.system('mv 2ML* KaKs/{}' .format(gene))
		os.system('mv 2NG* KaKs/{}' .format(gene))
		os.system('mv 4fold.nuc KaKs/{}'.format(gene))

def make_dnds_csv(codeml_file, all_gene_dict, gene):
	results = codeml.read(codeml_file)
	done = []
	if not os.path.exists('{}_dnds.csv' .format(gene)):
		dnds_file = open('{}_dnds.csv' .format(gene),'a+')
		dnds_file.write('dN/dS ratios pairwise,' + ','.join([x[1] for x in sorted(all_gene_dict.keys()) if x[4] == gene]))
		dnds_file.close()
	for tupl in sorted(all_gene_dict.keys()):
		dnds_file = open('{}_dnds.csv' .format(gene),'a+')
		if tupl[4] == gene:
			dnds_file.write('\n' + tupl[1] + ',')
			done.append(tupl[1])
			for other_tpl in sorted(all_gene_dict.keys()):
				if other_tpl[1] != tupl[1] and other_tpl[4] == gene and other_tpl[1] in done:
					dnds_file.write(str(results['pairwise'][tupl[1]][other_tpl[1]]['omega']) + ',')
				#elif other_tpl[1] == tupl[1] and other_tpl[4] == gene:
				#	dnds_file.write(',\n')
				else: continue
		dnds_file.close()

def write_evol_characteristics(file = '', gene = '', cai = 'NA', rscu = {}):
	if not os.path.exists('Codon_Analyses.csv'):
		final_file = open('Codon_Analyses.csv', 'a+')
		final_file.write('Names,Accession,Taxonomy,State,Gene,CAI,' + ','.join([x for x in sorted(rscu.keys())]) + '\n')
		final_file.close()
	final_file = open('Codon_Analyses.csv', 'a+')
	final_file.write(file[0] + ',' + file[1] + ',' + file[2] + ',' + file[3] + ',' + file[4] + ',' + str(cai) + ',' + ','.join([str(rscu[x]) for x in sorted(rscu.keys())]) + '\n')
	final_file.close()

def write_output(name = '', acc_number = '', family = '', mito_char = '', gene_char = '', state = '', assembly = '', sequencing = ''):
	if not os.path.exists('Final_Results.csv'):
		final_file = open('Final_Results.csv', 'a+')
		final_file.write('Name,Accession,Content,State,Assembly_Method,Technology,Taxonomy,Length,GC_content,%G,%C,%A,%T,GC_skew,AT_skew,Start,Stop,%G1,%C1,%A1,%T1,%G2,%C2,%A2,%T2,%G3,%C3,%A3,%T3\n')
		final_file.write(name + ',' + acc_number + ',' + 'Mitogenome' + ',' + state + ',' + assembly + ',' + sequencing + ',' + family + ',' + str(mito_char[0]) + ',' + str(mito_char[1]) + ',' + str(mito_char[2]) + ',' + str(mito_char[3]) + ',' + str(mito_char[4]) + ',' + str(mito_char[5]) + ',' + str(mito_char[6]) + ',' + str(mito_char[7]) + '\n')
		for gene in gene_char:
			final_file.write(name + ',' + acc_number + ','  +  gene  + ',' + state + ',' + assembly + ',' + sequencing + ',' + family + ',' + str(gene_char[gene][0]) + ',' + str(gene_char[gene][1]) + ',' + str(gene_char[gene][2]) + ',' + str(gene_char[gene][3]) + ',' + str(gene_char[gene][4]) + ',' + str(gene_char[gene][5]) + ',' + str(gene_char[gene][6]) + ',' + str(gene_char[gene][7]) + ',' + str(gene_char[gene][8]) + ',' + str(gene_char[gene][9]) + ',' + str(gene_char[gene][10]) + ',' + str(gene_char[gene][11]) + ',' + str(gene_char[gene][12]) + ',' + str(gene_char[gene][13]) + ',' + str(gene_char[gene][14]) + ',' + str(gene_char[gene][15]) + ',' + str(gene_char[gene][16]) + ',' + str(gene_char[gene][17]) + ',' + str(gene_char[gene][18]) + ',' + str(gene_char[gene][19]) + ',' + str(gene_char[gene][20]) + ',' + str(gene_char[gene][21]) + '\n')
		final_file.close()
	else:
		final_file = open('Final_Results.csv', 'a+')
		final_file.write(name + ',' + acc_number + ',' + 'Mitogenome'  + ',' + state + ',' + assembly + ',' + sequencing +  ',' + family + ',' + str(mito_char[0]) + ',' + str(mito_char[1]) + ',' + str(mito_char[2]) + ',' + str(mito_char[3]) + ',' + str(mito_char[4]) + ',' + str(mito_char[5]) + ',' + str(mito_char[6]) + ',' + str(mito_char[7]) + '\n')
		if gene_char != []:
			for gene in gene_char:
				final_file.write(name + ',' + acc_number + ','  +  gene + ',' + state + ',' + assembly + ',' + sequencing  + ',' + family + ',' + str(gene_char[gene][0]) + ',' + str(gene_char[gene][1]) + ',' + str(gene_char[gene][2]) + ',' + str(gene_char[gene][3]) + ',' + str(gene_char[gene][4]) + ',' + str(gene_char[gene][5]) + ',' + str(gene_char[gene][6]) + ',' + str(gene_char[gene][7]) + ',' + str(gene_char[gene][8]) + ',' + str(gene_char[gene][9]) + ',' + str(gene_char[gene][10]) + ',' + str(gene_char[gene][11]) + ',' + str(gene_char[gene][12]) + ',' + str(gene_char[gene][13]) + ',' + str(gene_char[gene][14]) + ',' + str(gene_char[gene][15]) + ',' + str(gene_char[gene][16]) + ',' + str(gene_char[gene][17]) + ',' + str(gene_char[gene][18]) + ',' + str(gene_char[gene][19]) + ',' + str(gene_char[gene][20]) + ',' + str(gene_char[gene][21]) + '\n')
		final_file.close()

def write_log(unverified = '', missing_gene = '', no_features = ''):
	log_file = open('log','a+')
	if no_features != []:
		log_file.write('###  THE FOLLOWING SEQUENCES DO NOT HAVE FEATURES:  ###\n\n')
		log_file.write('\n'.join(map(str,no_features)) + '\n')
	if unverified != []:
		log_file.write('###  THE FOLLOWING SEQUENCES ARE UNVERIFIED:  ###\n\n')
		log_file.write('\n'.join(map(str,unverified)) + '\n')
	if missing_gene != []:
		log_file.write('###  THE FOLLOWING GENES ARE MISSING:  ###\n\n')
		log_file.write('\n'.join(map(str,missing_gene)) + '\n')
	log_file.close()
	print('\nAll done! :)')

if __name__ == "__main__":
	args = parse_args()
	main(args)