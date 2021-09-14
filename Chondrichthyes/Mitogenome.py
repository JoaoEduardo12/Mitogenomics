from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class Mitogenome:
    def __init__(self, record, sequence, label, name, gen_code):
        '''Class construction with instance objects of the class Mitogenome'''
        self.record = record # the file name used
        self.sequence = sequence # complete mitogenome sequence
        self.label = label # all gene sequences
        self.name = name # a list with all taxonomic classifications
        self.gen_code = gen_code

    def get_id(self):
        return self.record

    def get_sequence(self):
        if self.label == "CDS":
            if (len(self.sequence) + 1) % 3 == 0:
                self.sequence += 'A'
            elif (len(self.sequence) + 2) % 3 == 0:
                self.sequence += 'AA'
        return self.sequence.upper()

    def get_aminoacid_sequence(self):
        aminoacid_sequence = self.sequence.translate(table = self.gen_code)
        aminoacid_sequence = str(aminoacid_sequence)
        if aminoacid_sequence[-1] == "*":
            aminoacid_sequence = aminoacid_sequence[:-1]
        return aminoacid_sequence

    def get_length(self):
        return len(self.sequence)

    def get_gc_content(self):
        return GC(self.sequence)

    def get_nucleotide_frequencies(self):
        nucleotide_freq = [str(self.sequence).upper().count(x)/ len(self.sequence) * 100 for x in ["A","T","G","C"]]
        return nucleotide_freq

    def get_skews(self):
        nucleotide_freq = self.get_nucleotide_frequencies()
        if (nucleotide_freq[0] + nucleotide_freq[1]) != 0 and (nucleotide_freq[2] + nucleotide_freq[3]) != 0:
            at_sk = (nucleotide_freq[0] - nucleotide_freq[1]) / (nucleotide_freq[0] + nucleotide_freq[1])
            gc_sk = (nucleotide_freq[2] - nucleotide_freq[3]) / (nucleotide_freq[2] + nucleotide_freq[3])
        else: at_sk, gc_sk = 0
        skew_vector = [at_sk, gc_sk]
        return skew_vector

    def get_cds_freq(self):
        if self.label == "CDS":
            codon_list = [str(self.sequence)[base:base+3] for base in range(0, len(str(self.sequence)), 3)]
            cds_freq = {}
            for i in range(3):
                cds_freq[i] = [''.join([base[i] for base in codon_list if len(base) == 3]).upper().count(x)/len(codon_list)*100 for x in ["A","T","G","C"]]
            return cds_freq
        else: return {0: ["","","",""], 1: ["","","",""], 2: ["","","",""]}

    def get_start_and_stop(self):
        if self.label == "CDS":
            return [self.sequence[0:3], self.sequence[-3::]]
        else: return ["",""]

    def get_molecular_weight(self):
        aminoacid_sequence = self.get_aminoacid_sequence()
        if self.label == "CDS" and "X" not in aminoacid_sequence and "B" not in aminoacid_sequence:
            return ProteinAnalysis(str(aminoacid_sequence)).molecular_weight()
        else: return ""

    def get_aromaticity(self):
        aminoacid_sequence = self.get_aminoacid_sequence()
        if self.label == "CDS" and "X" not in aminoacid_sequence and "B" not in aminoacid_sequence:
            return ProteinAnalysis(str(aminoacid_sequence)).aromaticity()
        else: return ""

    def get_instability_index(self):
        aminoacid_sequence = self.get_aminoacid_sequence()
        if self.label == "CDS" and "X" not in aminoacid_sequence and "B" not in aminoacid_sequence:
            return ProteinAnalysis(str(aminoacid_sequence)).instability_index()
        else: return ""

    def get_isoelectric_point(self):
        aminoacid_sequence = self.get_aminoacid_sequence()
        if self.label == "CDS" and "X" not in aminoacid_sequence and "B" not in aminoacid_sequence:
            return ProteinAnalysis(str(aminoacid_sequence)).isoelectric_point()
        else: return ""

    def get_gravy(self):
        aminoacid_sequence = self.get_aminoacid_sequence()
        if self.label == "CDS" and "X" not in aminoacid_sequence and "B" not in aminoacid_sequence:
            return ProteinAnalysis(str(aminoacid_sequence)).gravy()
        else: return ""
