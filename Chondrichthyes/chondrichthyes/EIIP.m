function[seq] = EIIP(seq)
% Input: Sequence of nucleotides A, T, G, C
% Output: Sequence in numerical representation electron-ion interaction
% potential.
% ---------------------------------------------------------------------
% seq = EIIP(seq)

n=find(seq == 65);   %ASCII code for A
seq(n)=0.1260;

n=find(seq == 67);   %ASCII code for C
seq(n)=0.1340;

n=find(seq == 71);   %ASCII code for G
seq(n)=0.0806;

n=find(seq == 84);   %ASCII code for T
seq(n)=0.1335;

end