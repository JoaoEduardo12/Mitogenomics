clc
clear all

%seq1 = fastaread('A3IRAN/A3IRAN_cds.fasta');
%seq2 = fastaread('BIV5907/BIV5907_cds.fasta');
seq1 = fastaread('BIV5002/BIV5002_cds.fasta');
seq2 = fastaread('BIV5255/BIV5255_cds.fasta');
%seq1 = fastaread('BIV5493/BIV5493_cds.fasta');
%seq2 = fastaread('BIV6668/BIV6668_cds.fasta');

% Gene concatenation
seq1_H = [seq1(1).Sequence seq1(2).Sequence seq1(3).Sequence seq1(4).Sequence seq1(5).Sequence seq1(6).Sequence seq1(9).Sequence seq1(12).Sequence seq1(13).Sequence];
seq1_L = [seq1(7).Sequence seq1(8).Sequence seq1(10).Sequence seq1(11).Sequence];
seq2_H = [seq2(1).Sequence seq2(2).Sequence seq2(3).Sequence seq2(4).Sequence seq2(5).Sequence seq2(6).Sequence seq2(9).Sequence seq2(12).Sequence seq2(13).Sequence];
seq2_L = [seq2(7).Sequence seq2(8).Sequence seq2(10).Sequence seq2(11).Sequence];

% Calculate synonymous and nonsynonymous substitution rates for the H-strand
[score,alignment] = nwalign(nt2aa(seq1_H,'GeneticCode', 5),nt2aa(seq2_H,'GeneticCode', 5));
alg1 = seqinsertgaps(seq1_H,alignment(1,:));
alg2 = seqinsertgaps(seq2_H,alignment(3,:));

[dn_H,ds_H] = dnds(alg1,alg2);
KaKs_H = dn_H/ds_H;

% Calculate synonymous and nonsynonymous substitution rates for the L-strand
[score,alignment] = nwalign(nt2aa(seq1_L,'GeneticCode', 5),nt2aa(seq2_L,'GeneticCode', 5));
alg1 = seqinsertgaps(seq1_L,alignment(1,:));
alg2 = seqinsertgaps(seq2_L,alignment(3,:));

[dn_L,ds_L] = dnds(alg1,alg2);
KaKs_L = dn_L/ds_L;

% Plot

figure('Position', [280 120 1500 800])
bar([KaKs_H(:), KaKs_L(:)], 'grouped');
ylabel('Ka / Ks')
xlabel('L vs H strand')
title('Ka/Ks');
ax = gca;
ax.XTickLabel = {'Ka/Ks H strand','Ka/Ks L strand'};

figure('Position', [280 120 1500 800])
bar([dn_H(:), ds_H(:), dn_L(:), ds_L(:)], 'grouped');
ylabel('Ka / Ks')
xlabel('L vs H strand')
title('Ka/Ks');
ax = gca;
ax.XTickLabel = {'Ka H strand','Ks H strand','Ka L strand', 'Ks L strand'};

