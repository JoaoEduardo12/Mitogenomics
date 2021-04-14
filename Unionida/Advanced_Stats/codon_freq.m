clc
clear all

atp6 = fastaread('atp6.fas');
atp8 = fastaread('atp8.fas');
cob = fastaread('cob.fas');
cox1 = fastaread('cox1.fas');
cox2 = fastaread('cox2.fas');
cox3 = fastaread('cox3.fas');
nad1 = fastaread('nad1.fas');
nad2 = fastaread('nad2.fas');
nad3 = fastaread('nad3.fas');
nad4 = fastaread('nad4.fas');
nad4l = fastaread('nad4l.fas');
nad5 = fastaread('nad5.fas');
nad6 = fastaread('nad6.fas');

final_seq = '';
for i = 1:121
    final_seq = [final_seq atp6(i).Sequence atp8(i).Sequence cob(i).Sequence cox1(i).Sequence cox2(i).Sequence cox3(i).Sequence nad1(i).Sequence nad2(i).Sequence nad3(i).Sequence nad4(i).Sequence nad4l(i).Sequence nad5(i).Sequence nad6(i).Sequence];
end

figure(1)
count = my_codoncount(final_seq,'figure',true,'geneticcode','Invertebrate Mitochondrial');
title('Codon Frequency')