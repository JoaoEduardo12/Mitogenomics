clc
clear all
% codon usage
% dn/ds ratio

%file_list = dir();
%file_list = struct2cell(file_list);

%for i=1:length(dir)-2
 %   cd (file_list{1,i+2})
  %  if isfile('*_cds.fasta')
   %     biv = file_list{1,i};
      %  biv = fastaread('*_cds.fasta');
       % disp('YAAAAS')
   %else
       % disp('NOOO')
   % end
   % cd ..
%end

A3IRAN = fastaread('A3IRAN/A3IRAN_cds.fasta');
BIV5907 = fastaread('BIV5907/BIV5907_cds.fasta');
BIV6795 = fastaread('BIV6795/BIV6795_cds.fasta');
BIV5481 = fastaread('BIV5481/BIV5481_cds.fasta');
BIV5523 = fastaread('BIV5523/BIV5523_cds.fasta');
BIV5002 = fastaread('BIV5002/BIV5002_cds.fasta');
BIV6668 = fastaread('BIV6668/BIV6668_cds.fasta');
BIV5255 = fastaread('BIV5255/BIV5255_cds.fasta');
BIV5493 = fastaread('BIV5493/BIV5493_cds.fasta');
BIV5536 = fastaread('BIV5536/BIV5536_cds.fasta');
BIV5573 = fastaread('BIV5573/BIV5573_cds.fasta');
BIV5578 = fastaread('BIV5578/BIV5578_cds.fasta');
BIV5565 = fastaread('BIV5565/BIV5565_cds.fasta');

final_seq = '';
for i = 1:13
    final_seq = [final_seq A3IRAN(i).Sequence BIV5907(i).Sequence BIV6795(i).Sequence BIV5481(i).Sequence BIV5523(i).Sequence BIV5002(i).Sequence BIV6668(i).Sequence BIV5255(i).Sequence BIV5493(i).Sequence BIV5536(i).Sequence BIV5573(i).Sequence BIV5578(i).Sequence];
end

figure(1)
count = my_codoncount(final_seq,'figure',true,'geneticcode','Invertebrate Mitochondrial');
%y(1,:) = cell2mat(struct2cell(count));
%count = struct2table(count);
%x(1,:) = count.Properties.VariableNames;
%heatmapHandle = heatmap(count,x,y,'ColorMap', jet(100));
title('Codon Frequency')


%%%%%%%%%%%%%%% Calculating base

KaKs = zeros(1,numel(A3IRAN));
for iCDS = 1:numel(A3IRAN)
        % align aa sequences of corresponding genes
        [score,alignment] = nwalign(nt2aa(A3IRAN(iCDS).Sequence,'GeneticCode', 5),nt2aa(BIV5907(iCDS).Sequence,'GeneticCode', 5));
        seq1 = seqinsertgaps(A3IRAN(iCDS).Sequence,alignment(1,:));
        seq2 = seqinsertgaps(BIV5907(iCDS).Sequence,alignment(3,:));

        % Calculate synonymous and nonsynonymous substitution rates
        [dn,ds] = dnds(seq1,seq2);
        KaKs(iCDS) = dn/ds;
end

KaKs_2 = zeros(1,numel(A3IRAN));
for iCDS = 1:numel(A3IRAN)
        % align aa sequences of corresponding genes
        [score,alignment] = nwalign(nt2aa(BIV6668(iCDS).Sequence,'GeneticCode', 5),nt2aa(BIV5255(iCDS).Sequence,'GeneticCode', 5));
        seq1 = seqinsertgaps(BIV6668(iCDS).Sequence,alignment(1,:));
        seq2 = seqinsertgaps(BIV5255(iCDS).Sequence,alignment(3,:));

        % Calculate synonymous and nonsynonymous substitution rates
        [dn,ds] = dnds(seq1,seq2);
        KaKs_2(iCDS) = dn/ds;
end

KaKs_3 = zeros(1,numel(A3IRAN));
for iCDS = 1:numel(A3IRAN)
        % align aa sequences of corresponding genes
        [score,alignment] = nwalign(nt2aa(BIV5493(iCDS).Sequence,'GeneticCode', 5),nt2aa(BIV5536(iCDS).Sequence,'GeneticCode', 5));
        seq1 = seqinsertgaps(BIV5493(iCDS).Sequence,alignment(1,:));
        seq2 = seqinsertgaps(BIV5536(iCDS).Sequence,alignment(3,:));

        % Calculate synonymous and nonsynonymous substitution rates
        [dn,ds] = dnds(seq1,seq2);
        KaKs_3(iCDS) = dn/ds;
end

KaKs_4 = zeros(1,numel(A3IRAN));
for iCDS = 1:numel(A3IRAN)
        % align aa sequences of corresponding genes
        [score,alignment] = nwalign(nt2aa(BIV5573(iCDS).Sequence,'GeneticCode', 5),nt2aa(BIV5578(iCDS).Sequence,'GeneticCode', 5));
        seq1 = seqinsertgaps(BIV5573(iCDS).Sequence,alignment(1,:));
        seq2 = seqinsertgaps(BIV5578(iCDS).Sequence,alignment(3,:));

        % Calculate synonymous and nonsynonymous substitution rates
        [dn,ds] = dnds(seq1,seq2);
        KaKs_4(iCDS) = dn/ds;
end

figure('Position', [280 120 1500 800])
bar(KaKs);
ylabel('Ka / Ks')
xlabel('genes')
title('Ka/Ks');
ax = gca;
ax.XTickLabel = {A3IRAN.Header};
%BIV5493, BIV5536
% plot Ka/Ks ratio for each gene
figure('Position', [280 120 1500 800])
bar([KaKs(:),KaKs_2(:),KaKs_3(:),KaKs_4(:)],'grouped');
ylabel('Ka / Ks')
xlabel('genes')
legend('Unionidae', 'Mulleridae','Margaritiferidae','Iridinidae');
title('Ka/Ks');
ax = gca;
ax.XTickLabel = {A3IRAN.Header};
% plot dotted line at threshold 1
%hold on
%line([0 numel(KaKs)+1],[1 1],'LineStyle', ':');
KaKs

%%%%%% Calculating sliding window

% ORF number corresponding to gene ENV

% align the two ENV genes

[cox1_1, cox1_2, dn1, ds1] = get_alignments(A3IRAN,BIV5907, 1);
[cox3_1, cox3_2, dn2, ds2] = get_alignments(A3IRAN,BIV5907, 2);
[atp6_1, atp6_2, dn3, ds3] = get_alignments(A3IRAN,BIV5907, 3);
[atp8_1, atp8_2, dn4, ds4] = get_alignments(A3IRAN,BIV5907, 4);
[nad4l_1, nad4l_2, dn5, ds5] = get_alignments(A3IRAN,BIV5907, 5);
[nad4_1, nad4_2, dn6, ds6] = get_alignments(A3IRAN,BIV5907, 6);
[nad6_1, nad6_2, dn7, ds7] = get_alignments(A3IRAN,BIV5907, 7);
[nad1_1, nad1_2, dn8, ds8] = get_alignments(A3IRAN,BIV5907, 8);
[nad5_1, nad5_2, dn9, ds9] = get_alignments(A3IRAN,BIV5907, 9);
[cob_1, cob_2, dn10, ds10] = get_alignments(A3IRAN,BIV5907, 10);
[nad2_1, nad2_2, dn11, ds11] = get_alignments(A3IRAN,BIV5907, 11);
[nad3_1, nad3_2, dn12, ds12] = get_alignments(A3IRAN,BIV5907, 12);
[cox2_1, cox2_2, dn13, ds13] = get_alignments(A3IRAN,BIV5907, 13);

% plot the Ka/Ks trends for the different window sizes
figure('Position', [20 20 2000 1200])
hold on
subplot(7,2,1)
plot(dn1./ds1, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('COX1')

subplot(7,2,2)
plot(dn2./ds2, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('COX3')

subplot(7,2,3)
plot(dn3./ds3, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('ATP6')

subplot(7,2,4)
plot(dn4./ds4, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('ATP8')

subplot(7,2,5)
plot(dn5./ds5, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('NAD4L')

subplot(7,2,6)
plot(dn6./ds6, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('NAD4')

subplot(7,2,7)
plot(dn7./ds7, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('NAD6')

subplot(7,2,8)
plot(dn8./ds8, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('NAD1')

subplot(7,2,9)
plot(dn9./ds9, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('NAD5')

subplot(7,2,10)
plot(dn10./ds10, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('COB')

subplot(7,2,11)
plot(dn11./ds11, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('NAD2')

subplot(7,2,12)
plot(dn12./ds12, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('NAD3')

subplot(7,2,13)
plot(dn13./ds13, 'b');
ylabel('Ka / Ks')
ylim([0 10])
xlabel('sliding window (starting codon)') % ORF number corresponding to gene ENV
title('COX2')

