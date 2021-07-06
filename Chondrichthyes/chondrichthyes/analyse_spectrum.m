clc
clear all

atp6 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/atp6.fn');
atp8 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/atp8.fn');
cob = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/cob.fn');
cox1 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/cox1.fn');
cox2 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/cox2.fn');
cox3 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/cox3.fn');
nad1 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/nad1.fn');
nad2 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/nad2.fn');
nad3 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/nad3.fn');
nad4 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/nad4.fn');
nad4l = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/nad4l.fn');
nad5 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/nad5.fn');
nad6 = fastaread('/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Corrected_Genes/nad6.fn');


atp6_seq = EIIP(atp6(1).Sequence);
