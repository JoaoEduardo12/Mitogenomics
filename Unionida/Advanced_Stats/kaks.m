
clear all

gene = fastaread('nad6.fas');

ka_ks = zeros(120,120);
ka = zeros(120,120);
ks = zeros(120,120);

for i = 1:119
    for j = i:119
        if strcmp(gene(i).Header,gene(j).Header)
            continue
        else
            [score,alignment] = nwalign(nt2aa(gene(i).Sequence,'GeneticCode', 5, 'ACGTOnly',false),nt2aa(gene(j).Sequence,'GeneticCode', 5, 'ACGTOnly',false));
            seq1 = seqinsertgaps(gene(i).Sequence,alignment(1,:));
            seq2 = seqinsertgaps(gene(j).Sequence,alignment(3,:));
            [dn,ds] = dnds(seq1,seq2,'method','PBL');
            ka(i,j) = dn;
            ks(i,j) = ds;
            ka_ks(i,j) = dn/ds;
        end
    end
end

writetable(table(ka),'nad6_ka.xlsx')
writetable(table(ka_ks),'nad6_kaks.xlsx')
writetable(table(ks),'nad6_ks.xlsx')