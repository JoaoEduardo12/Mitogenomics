import sys, os
import argparse
import re
from Bio import SeqIO
from Bio import Entrez
from CAI import RSCU
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Phylo.PAML import codeml
from CAI.CAI import CAI
from Mitogenome import Mitogenome
Entrez.email='joao.edu.t@hotmail.com'
api_key = '105bc09b7b34973aa305f31be73aecbc4109'

def parse_args():
    """Parsing command line arguments"""
    parser = argparse.ArgumentParser(prog="deep_mitogenomics")
    parser = argparse.ArgumentParser(description='Given a search term in NCBI, reads all available sequences and extracts multiple statistics')
    parser.add_argument('-search', '--Search', help = 'Give this command a file with the search term in NCBI (right side of the page in "Search Details, when you query an entry")')
    parser.add_argument('-out', '--OutFile', help = "Give this command the output file name")
    parser.add_argument('-common_taxonomy', '--CommTax', help = 'Give this command the highest recent common taxonomy classification within search (optional)')
    parser.add_argument('-target_taxonomy', '--TargetTax', help = 'Give this command the taxonomy in which you wish to make distinctions (i.e. will add an extra column dictating this taxonomy rank. Choose between family and subfamily')
    parser.add_argument('-genetic_code', '--GeneCode', help = 'Give this command the genetic code of the input sequences. Defaults to 11', default = 11)
    parser.add_argument('-gene_download', "--DownloadGeneDir", help = '(Optional) Give this command a directory name where the mitochondrial gene files will be saved', default = False)
    parser.add_argument('-outgroup', "--Outgroup", help = '(Optional only if gene_download is selected) Give this command an outgroup species fasta file with gene features from NCBI', default = False)
    parser.add_argument('-mito_download', "--DownloadMitoDir", help = '(Optional) Give this command a directory name where the complete mitogenomes will be saved', default = False)
    parser.add_argument('-clean', "--Clean", help = '(Optional) Give this command the argument to remove sequences which have at least one gene with gaps', default = True)
    parser.add_argument('-dnds', "--DnDs", help = '(Optional) if True will calculate dnds values of each protein coding gene, if the command --codeml_model is not given, the M0 model (the most simple) will be computed', default = False)
    parser.add_argument('-model', "--Model", help = '(Optional if dnds = True, otherwise it can be ignored) Select one of the models available in codeml (consult its documentation), options include M0, M1, M2, M3, M4, M7, and M8. Note that some of the options require a phylogenetic tree', default = False)
    parser.add_argument('-codon_properties', "--CodonProperties", help = 'if codon_properties = True it will calculate codon properties, e.g., effective number of codon, codon adaptation index, and relative synonymous codon usage (default = False)', default = False)
    parser.add_argument('-cai',"--CAI", help="Give this command a highly expressed mitochondrial gene so that the CAI index may be calculated (beware that this will greatly increase the computation time)")
    parser.add_argument('-unique_species',"--UniqueSpecies", help="Give this command the option to keep only mitogenomes of one species at a time, so multiple species will not appear")
    parser.add_argument('-remove_accessions', "--Accessions",  help='Give this command a file with accession numbers marked with hashtag beforehand (anywhere)')
    args = parser.parse_args()
    return args

def main(args):
    """Uses the NCBI search utilities to retrieve sequences and genes inside"""
    if args.DownloadGeneDir:
        os.mkdir(args.DownloadGeneDir)
    if args.DownloadMitoDir:
        os.mkdir(args.DownloadMitoDir)
    unverified, no_features, no_genes = [], [], [] # list containing warnings throughout this script
    print('\nSearching in GenBank\n')
    handle = Entrez.esearch(db='nucleotide', retmax=100000, term= args.Search, idtype="acc") # Search term to retrieve sequences
    record = Entrez.read(handle)
    print('\n{} results found\n' .format(str(len(record['IdList']))))
    handle.close()
    print('\nExtracting features\n')
    i = 0
    cds_dict = {}
    species_list = []
    for rec in record['IdList']:
        i += 1
        print(' ' + str(round(i/len(record['IdList'])*100, 2)) + '%', end = "\r")
        unverified, no_features, no_genes, cds_dict = extract_features(rec, unverified, no_features, no_genes, cds_dict, species_list, args) # This is the main function that will extract complete mitogenomes and its genes, in addition of performing analyses on them
    if args.Outgroup:
        process_outgroup(args.Outgroup, args)
    if args.CodonProperties:
        calculate_codon_properties(cds_dict, args) # These functions from here on, will need all the gathered information from all records, as it implies comparative analyses
    #if args.DnDs:
    #    calculate_dnds(Mitogenome, args.Model)

def extract_features(record, unverified, no_features, no_genes, cds_dict, species_list, args):
    """
    This function will gather all the information into the Mitogenome class, and process
    other information as well, as attributed in the user options.
    """
    if args.Accessions: # This is to "manually" filter accession numbers that the user does not want
        records_to_exclude = filter_accessions(args.Accessions)
    else: records_to_exclude = []
    mitogenome_characteristics = []
    gene_characteristics = {}
    gene_list = [] # this is an empty list that will have all the gene names, to later confirm that a mitogenome has all of its genes
    with Entrez.efetch(db="nucleotide", rettype = "gb", retmode = "text", id = record) as handle:
        if record not in records_to_exclude:
            # the following 2 lines of code is to read the complete mitogenome
            full_sequence = Entrez.efetch(db = "nucleotide", id = record, rettype = 'fasta', retmode = 'text')
            full_sequence = full_sequence.read()
            if args.Clean:
                clean = filter_sequence(full_sequence, record)
            if clean:
                label = "Mitogenomes" # this is just for the class Mitogenome to know that it is in fact a complete mitogenome
                name = "Mitogenomes"
                # using the class Mitogenome to store the sequence and calculate its genomic statistics
                mitogenome_characteristics = Mitogenome(record, full_sequence, label, name, args.GeneCode)
                if args.DownloadMitoDir:
                    write_sequence_files(name, record, mitogenome_characteristics, args.DownloadMitoDir, False, True, args.GeneCode)
                for rec in SeqIO.parse(handle, "gb"): # This now will take care of the information about taxonomy, among others, and of course of the genes inside the record
                    mitogenome_switch = True
                    if args.UniqueSpecies:
                        if rec.annotations["organism"] in species_list:
                            mitogenome_switch = False
                        else: 
                            species_list.append(rec.annotations["organism"])
                    if mitogenome_switch == True:
                        taxonomy_names = ";".join(rec.annotations["taxonomy"])
                        # These lines of codes take the info directly from NCBI, whether they are complete, incomplete or unverified
                        if "UNVERIFIED" in rec.description:
                            unverified.append(rec.description)
                            state = "unverified"
                        elif "partial" in rec.description:
                            state = "partial"
                        else:
                            state = "complete" # it will later create a column in the csv with this information
                        name_list = []
                        if rec.features:
                            for feature in rec.features:
                                if feature.type == "CDS": # if the information is a protein coding gene
                                    label = "CDS" #this will tell the created class Mitogenome, that the sequence that is being calculated upon, is this
                                    if "gene" in feature.qualifiers.keys():
                                        name = feature.qualifiers["gene"][0].strip('') # gene name
                                    else:
                                        name = feature.qualifiers["product"][0].strip('') 
                                    gene_sequence = feature.location.extract(rec).seq # gene sequence
                                    name = gene_name_correction(name, gene_characteristics)
                                    gene_characteristics[name] = Mitogenome(record, gene_sequence, label, name, args.GeneCode)
                                    # a dictionary will ber created with tuple values denoting every gene from every sequence, we want to store this information to calculate metrics that require comparison later on
                                    cds_dict[(rec.annotations['organism'],rec.id,taxonomy_names.split(args.CommTax + ";")[1],state,name)] = gene_characteristics[name].get_sequence() # we are telling the class to return the sequence, this is needed because inside the class, further processing is done in order to clean the sequence
                                elif feature.type == "tRNA": # transfer RNA
                                    # Same procedures through and through
                                    label = "tRNA"
                                    gene_sequence = feature.location.extract(rec).seq
                                    if "product" in feature.qualifiers.keys():
                                        name = feature.qualifiers['product'][0].strip('')
                                        name = gene_name_correction(name, gene_characteristics)
                                        gene_characteristics[name] = Mitogenome(record, gene_sequence, label, name, args.GeneCode)
                                elif feature.type == "rRNA":      # ditto
                                    label = "rRNA"
                                    gene_sequence = feature.location.extract(rec).seq
                                    if "product" in feature.qualifiers.keys():
                                        name = feature.qualifiers['product'][0].strip('')
                                        name = gene_name_correction(name, gene_characteristics)
                                        gene_characteristics[name] = Mitogenome(record, gene_sequence, label, name, args.GeneCode)
                                elif feature.type == "source": # some features have the label source, we want to ignore this, so the instance gene_characteristics will be initialized to a empty list
                                    gene_characteristics = {}
                                if args.DownloadGeneDir and gene_characteristics != {}:
                                    if name not in name_list:
                                        write_sequence_files(name, rec.id, gene_characteristics, args.DownloadGeneDir, True, True, args.GeneCode)
                                    name_list.append(name)
                        else:
                            no_features.append(rec.description) # if this line of code runs, it means that at least in the genbank file there are no information regarding genes
                            gene_characteristics = {}
                        # we want also to obtain information about the sequencing method (sanger, ngs, tgs) and the assembly programs used, at least the ones who do have that information
                        assembly_method = ""
                        sequencing_technology = ""
                        if "structured_comment" in rec.annotations:
                            if "Assembly-Data" in rec.annotations["structured_comment"].keys():
                                if "Assembly Method" in rec.annotations['structured_comment']['Assembly-Data'].keys():
                                    assembly_method = rec.annotations['structured_comment']['Assembly-Data']['Assembly Method']
                                if 'Sequencing Technology' in rec.annotations['structured_comment']['Assembly-Data'].keys():
                                    sequencing_technology = rec.annotations['structured_comment']['Assembly-Data']['Sequencing Technology']
                        # In the following line of code the taxonomy_names variable will be cut to include the user defined common taxonomy of all sequences, leaving the remaining values stored for unique taxonomic groups
                        if args.CommTax:
                            taxonomy_names = taxonomy_names.split(args.CommTax + ";")[1]
                        # In this line, we are just proceeding carefully to know if the mitogenome has all genes available, if it does not the tag name state is changed to partial
                        if gene_characteristics != {}:
                            no_genes = gene_check(gene_list, no_genes, rec.annotations["organism"], taxonomy_names)
                        else: state = "partial"
                        # At last we will be writing all information gathered in the output file, along the entire script!
                        write_output(name = rec.annotations["organism"], acc_number = rec.id, taxonomy_group = taxonomy_names, mitogenome_char = mitogenome_characteristics, gene_char = gene_characteristics, state = state, assembly = assembly_method, sequencing = sequencing_technology, args= args)
    return unverified, no_features, no_genes, cds_dict

def process_outgroup(file, args):
    gene_list = []
    pcg = ['cox1','cox3','atp6','atp8','nad4l','nad4','nad6','nad1','nad5','cob','nad2','nad3','cox2']
    for rec in SeqIO.parse(file, "fasta"):
        accession = rec.id.split("|")[1]
        accession = accession.split("_gene")[0]
        gene_name = rec.description.split("gene=")[1]
        gene_name = gene_name.split("]")[0]
        gene_name = gene_name_correction(gene_name, gene_list)
        gene_list.append(gene_name)
        if gene_name in pcg:
            write_sequence_files(gene_name, accession, rec.seq, args.DownloadGeneDir, True, False, args.GeneCode)
        else: write_sequence_files(gene_name, accession, rec.seq, args.DownloadGeneDir, False, False, args.GeneCode)

def filter_accessions(accessions):
	records = []
	with open(accessions) as f:
		lines = f.readlines()
		f.close()
	for line in lines:
		if '#' in line:
			records.append(line.strip('\n')[1:])
	return records

def filter_sequence(seq, id):
	clean = True
	if seq.count('N') > 7:
		print('\n{} is no good, removing it from the dataset' .format(id))
		clean = False
	return clean

def gene_check(gene_list, no_genes, sample, tax):
	complete_list = ['atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6','rrnL','rrnS','tRNA-Phe','tRNA-Val','tRNA-Leu','tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Ala','tRNA-Asn','tRNA-Cys','tRNA-Tyr','tRNA-Ser','tRNA-Asp','tRNA-Lys','tRNA-Gly','tRNA-Arg','tRNA-His','tRNA-Ser','tRNA-Leu','tRNA-Glu','tRNA-Thr','tRNA-Pro','tRNA-Ser2','tRNA-Leu2']
	for gene in gene_list:
		if not gene in complete_list:
			gene_report = 'The gene {} is not present in {} (Taxonomy: {})' .format(gene,sample,tax)
			no_genes.append(gene_report)
		else: continue

def gene_name_correction(name, gene_characteristics):
    if re.fullmatch("(MT)?atp( )?(-)?(synth)?(ase)?( )?(subunit)?( )?(F0)?( )?(subunit)?( )?(-)?(6|vi)",name, re.IGNORECASE):
        new_name = "atp6"
    elif re.fullmatch("(MT)?atp( )?(-)?(synth)?(ase)?( )?(subunit)?( )?(F0)?( )?(subunit)?( )?(-)?(8|viii)", name, re.IGNORECASE):
        new_name = "atp8"
    elif re.fullmatch("[(]?c(ytochrome)?( )?(c)?( )?(oxidase)?( )?(subunit)?( )?(o)?(0)?( )?(x)?(-)?( )?(1|I)[)]?", name, re.IGNORECASE):
        new_name = "cox1"
    elif re.fullmatch("[(]?c(ytochrome)?( )?(c)?( )?(oxidase)?( )?(subunit)?( )?(o)?(0)?( )?(x)?(-)?( )?(2|II)[)]?", name, re.IGNORECASE):
        new_name = "cox2"
    elif re.fullmatch("[(]?c(ytochrome)?( )?(c)?( )?(oxidase)?( )?(subunit)?( )?(o)?(0)?( )?(x)?(-)?( )?(3|III)[)]?", name, re.IGNORECASE):
        new_name = "cox3"
    elif re.fullmatch("c(ob)?(yt)?(ochrome)?(_)?(-)?( )?(y)?( )?(t)?( )?(b)?", name, re.IGNORECASE):
        new_name = "cob"
    elif re.fullmatch("n(a)?(d)?( )?(-)?(h)?( )?(dehydrogenase)?( )?(subunit)?( )?(-)?(1|I)", name, re.IGNORECASE):
        new_name = "nad1"
    elif re.fullmatch("n(a)?(d)?( )?(-)?(h)?( )?(dehydrogenase)?( )?(subunit)?( )?(-)?(2|II)", name, re.IGNORECASE):
        new_name = "nad2"
    elif re.fullmatch("n(a)?(d)?( )?(-)?(h)?( )?(dehydrogenase)?( )?(subunit)?( )?(-)?(3|III)", name, re.IGNORECASE):
        new_name = "nad3"
    elif re.fullmatch("n(a)?(d)?( )?(-)?(h)?( )?(dehydrogenase)?( )?(subunit)?( )?(-)?(4|IV)", name, re.IGNORECASE):
        new_name = "nad4"
    elif re.fullmatch("n(a)?(d)?( )?(-)?(h)?( )?(dehydrogenase)?( )?(subunit)?( )?(-)?(4|IV)( )?(-)?l", name, re.IGNORECASE):
        new_name = "nad4l"
    elif re.fullmatch("n(a)?(d)?( )?(-)?(h)?( )?(dehydrogenase)?( )?(subunit)?( )?(-)?(5|V)", name, re.IGNORECASE):
        new_name = "nad5"
    elif re.fullmatch("n(a)?(d)?( )?(-)?(h)?( )?(dehydrogenase)?( )?(subunit)?( )?(-)?(6|VI)", name, re.IGNORECASE):
        new_name = "nad6"
    elif name in ['D-loop','D-loop','control region','d-loop','D loop','D Loop','D-Loop','D LOOP','D-LOOP','CONTROL REGION']:
        new_name = "D-loop"
    elif name in ['tRNA-Phe','tRNA-Val','tRNA-Leu','tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Ala','tRNA-Asn','tRNA-Cys','tRNA-Tyr','tRNA-Ser','tRNA-Asp','tRNA-Lys','tRNA-Gly','tRNA-Arg','tRNA-His','tRNA-Ser','tRNA-Leu','tRNA-Glu','tRNA-Thr','tRNA-Pro',"trnF","trnV","trnL1","trnI","trnQ","trnM","trnW","trnA","trnN","trnC","trnY","trnS1","trnD","trnK","trnG","trnR","trnH","trnS2","trnL2","trnT","trnE","trnP"]:
        new_name = name
    elif name in ['rrnL','16S','16S ribosomal RNA','16S Ribosomal RNA','l6S ribosomal RNA','l-rRNA','L-rRNA','l rRNA','L rRNA','l-rrna','l rrna','L RRNA','lrrna','l rRNA','lrRNA',"rrn16"]:
        new_name = "rrnL"
    elif name in ['rrnS','12S','12S ribosomal RNA','12S Ribosomal RNA','l2S ribosomal RNA','s-rRNA','S-rRNA','s rRNA','S rRNA','s-rrna','s rrna','S RRNA','srrna','s rRNA','srRNA',"rrn12"]:
        new_name = "rrnS"
    else:
        sys.exit("There seems to be a gene read that does not appear to be a typical mitochondrial gene: %s" % name)
    if new_name in gene_characteristics:
        new_name += "2"
    return new_name

def write_sequence_files(sequence_name, sequence_id, sequence, directory, gene_switch, nonoutgroup_switch, genetic_code):
    pcg = ['cox1','cox3','atp6','atp8','nad4l','nad4','nad6','nad1','nad5','cob','nad2','nad3','cox2']
    if nonoutgroup_switch:
        if gene_switch:
            with open(os.path.join(directory, sequence_name) + ".fas", "a+") as nuc_file:
                nuc_file.write(">%s\n%s\n" % (sequence_id, sequence[sequence_name].get_sequence()))
                nuc_file.close()
            if sequence_name in pcg:
                with open(os.path.join(directory, sequence_name) + ".faa", "a+") as aa_file:
                    aa_file.write(">%s\n%s\n" % (sequence_id, sequence[sequence_name].get_aminoacid_sequence()))
                    aa_file.close()
        else:
            with open(os.path.join(directory, sequence_name) + ".fas", "a+") as complete_file:
                complete_file.write(">%s\n%s\n" % (sequence_id, sequence.get_sequence()))
                complete_file.close()
    else:
        if gene_switch:
            with open(os.path.join(directory, sequence_name) + ".fas", "a+") as nuc_file:
                    nuc_file.write(">%s\n%s\n" % (sequence_id, sequence))
                    nuc_file.close()
            with open(os.path.join(directory, sequence_name) + ".faa", "a+") as aa_file:
                    aa_file.write(">%s\n%s\n" % (sequence_id, sequence.translate(table = genetic_code))) # meter aminoacidos nisto!
                    aa_file.close()
        else:
            with open(os.path.join(directory, sequence_name) + ".fas", "a+") as nuc_file:
                    nuc_file.write(">%s\n%s\n" % (sequence_id, sequence))
                    nuc_file.close()

def write_output(name = '', acc_number = '', taxonomy_group = "", mitogenome_char = '', gene_char = '', state = '', assembly = '', sequencing = '', args = ""):
    pcg = ['cox1','cox3','atp6','atp8','nad4l','nad4','nad6','nad1','nad5','cob','nad2','nad3','cox2']
    if not os.path.exists(args.OutFile + ".csv"):
        final_file = open(args.OutFile + '.csv', 'a+')
        final_file.write('Name,Accession,Content,State,Assembly_Method,Technology,Taxonomy,Length,GC_content,A,T,G,C,AT_skew,GC_skew,Start,Stop,A1,T1,G1,C1,A2,T2,G2,C2,G3,C3,A3,T3,Molecular_Weight,Aromaticity,Instability_Index,Isoelectric_Point,GRAVY\n')
        final_file.write(','.join([name, acc_number, "Mitogenome", state, assembly, sequencing, taxonomy_group, str(mitogenome_char.get_length()), str(mitogenome_char.get_gc_content()), str(mitogenome_char.get_nucleotide_frequencies()[0]), str(mitogenome_char.get_nucleotide_frequencies()[1]), str(mitogenome_char.get_nucleotide_frequencies()[2]), str(mitogenome_char.get_nucleotide_frequencies()[3]), str(mitogenome_char.get_skews()[0]), str(mitogenome_char.get_skews()[1]), "\n"]))
        for gene in gene_char:
            if gene in pcg:
                final_file.write(','.join([name, acc_number, gene, state, assembly, sequencing, taxonomy_group, str(gene_char[gene].get_length()), str(gene_char[gene].get_gc_content()), str(gene_char[gene].get_nucleotide_frequencies()[0]), str(gene_char[gene].get_nucleotide_frequencies()[1]), str(gene_char[gene].get_nucleotide_frequencies()[2]), str(gene_char[gene].get_nucleotide_frequencies()[3]), str(gene_char[gene].get_skews()[0]), str(gene_char[gene].get_skews()[1]), str(gene_char[gene].get_start_and_stop()[0]), str(gene_char[gene].get_start_and_stop()[1]), str(gene_char[gene].get_cds_freq()[0][0]), str(gene_char[gene].get_cds_freq()[0][1]), str(gene_char[gene].get_cds_freq()[0][2]), str(gene_char[gene].get_cds_freq()[0][3]), str(gene_char[gene].get_cds_freq()[1][0]), str(gene_char[gene].get_cds_freq()[1][1]), str(gene_char[gene].get_cds_freq()[1][2]), str(gene_char[gene].get_cds_freq()[1][3]), str(gene_char[gene].get_cds_freq()[2][0]), str(gene_char[gene].get_cds_freq()[2][1]), str(gene_char[gene].get_cds_freq()[2][2]), str(gene_char[gene].get_cds_freq()[2][3]), str(gene_char[gene].get_molecular_weight()), str(gene_char[gene].get_aromaticity()), str(gene_char[gene].get_instability_index()), str(gene_char[gene].get_isoelectric_point()), str(gene_char[gene].get_gravy()),"\n"]))
            else: final_file.write(','.join([name, acc_number, gene, state, assembly, sequencing, taxonomy_group, str(gene_char[gene].get_length()), str(gene_char[gene].get_gc_content()), str(gene_char[gene].get_nucleotide_frequencies()[0]), str(gene_char[gene].get_nucleotide_frequencies()[1]), str(gene_char[gene].get_nucleotide_frequencies()[2]), str(gene_char[gene].get_nucleotide_frequencies()[3]), str(gene_char[gene].get_skews()[0]), str(gene_char[gene].get_skews()[1]), "\n"]))
        final_file.close()
    else:
        final_file = open(args.OutFile + '.csv', 'a+')
        final_file.write(','.join([name, acc_number, "Mitogenome", state, assembly, sequencing, taxonomy_group, str(mitogenome_char.get_length()), str(mitogenome_char.get_gc_content()), str(mitogenome_char.get_nucleotide_frequencies()[0]), str(mitogenome_char.get_nucleotide_frequencies()[1]), str(mitogenome_char.get_nucleotide_frequencies()[2]), str(mitogenome_char.get_nucleotide_frequencies()[3]), str(mitogenome_char.get_skews()[0]), str(mitogenome_char.get_skews()[1]), "\n"]))
        for gene in gene_char:
            if gene in pcg:
                final_file.write(','.join([name, acc_number, gene, state, assembly, sequencing, taxonomy_group, str(gene_char[gene].get_length()), str(gene_char[gene].get_gc_content()), str(gene_char[gene].get_nucleotide_frequencies()[0]), str(gene_char[gene].get_nucleotide_frequencies()[1]), str(gene_char[gene].get_nucleotide_frequencies()[2]), str(gene_char[gene].get_nucleotide_frequencies()[3]), str(gene_char[gene].get_skews()[0]), str(gene_char[gene].get_skews()[1]), str(gene_char[gene].get_start_and_stop()[0]), str(gene_char[gene].get_start_and_stop()[1]), str(gene_char[gene].get_cds_freq()[0][0]), str(gene_char[gene].get_cds_freq()[0][1]), str(gene_char[gene].get_cds_freq()[0][2]), str(gene_char[gene].get_cds_freq()[0][3]), str(gene_char[gene].get_cds_freq()[1][0]), str(gene_char[gene].get_cds_freq()[1][1]), str(gene_char[gene].get_cds_freq()[1][2]), str(gene_char[gene].get_cds_freq()[1][3]), str(gene_char[gene].get_cds_freq()[2][0]), str(gene_char[gene].get_cds_freq()[2][1]), str(gene_char[gene].get_cds_freq()[2][2]), str(gene_char[gene].get_cds_freq()[2][3]), str(gene_char[gene].get_molecular_weight()), str(gene_char[gene].get_aromaticity()), str(gene_char[gene].get_instability_index()), str(gene_char[gene].get_isoelectric_point()), str(gene_char[gene].get_gravy()),"\n"]))
            else: final_file.write(','.join([name, acc_number, gene, state, assembly, sequencing, taxonomy_group, str(gene_char[gene].get_length()), str(gene_char[gene].get_gc_content()), str(gene_char[gene].get_nucleotide_frequencies()[0]), str(gene_char[gene].get_nucleotide_frequencies()[1]), str(gene_char[gene].get_nucleotide_frequencies()[2]), str(gene_char[gene].get_nucleotide_frequencies()[3]), str(gene_char[gene].get_skews()[0]), str(gene_char[gene].get_skews()[1]), "\n"]))
        final_file.close()

def calculate_codon_properties(cds_dict, args):
    print("\n\nPerforming codon analyses now\n")
    pcg = ['cox1','cox3','atp6','atp8','nad4l','nad4','nad6','nad1','nad5','cob','nad2','nad3','cox2']
    cds_dict = filter_ambiguous(cds_dict, pcg)
    reference = []
    i = 1
    cont = -1
    for tupl in cds_dict.keys():
        print(' ' + str(round(i/len(cds_dict.keys())*100, 2)) + '%', end = "\r")
        i += 1
        for gene in pcg:
            if gene in tupl:
                cai_value = "NA"
                sequence = cds_dict[tupl]
                sequence = sequence[:-3]
                rscu_list = []
                rscu_list.append(sequence)
                if args.CAI == True:
                    for other_tpl in cds_dict.keys():
                        cont += 1
                        if other_tpl != tupl:
                            if list(cds_dict.keys())[cont][4] == args.CAI:
                                reference.append(cds_dict[other_tpl])
                    cont = -1
                    cai_value = CAI(sequence, reference= reference, genetic_code= int(args.GeneCode))
                rscu_values = RSCU(rscu_list, genetic_code= int(args.GeneCode))
                write_codon_properties(file = tupl, cai = cai_value, rscu = rscu_values, args = args)

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

def write_codon_properties(file = "", cai = "", rscu = "", args = ""):
    if not os.path.exists(args.OutFile + "_codons.csv"):
        final_file = open(args.OutFile + "_codons.csv", "a+")
        final_file.write("Name,Accession,Taxonomy,State,Gene,CAI," + ",".join([x for x in sorted(rscu.keys())]) + "\n")
        final_file.close()
    final_file = open(args.OutFile + "_codons.csv", "a+")
    final_file.write(file[0] + ',' + file[1] + ',' + file[2] + ',' + file[3] + ',' + file[4] + ',' + str(cai) + ',' + ','.join([str(rscu[x]) for x in sorted(rscu.keys())]) + '\n')
    final_file.close()

if __name__ == "__main__":
	args = parse_args()
	main(args)
