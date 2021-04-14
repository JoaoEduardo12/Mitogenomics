function [x,y,dn,ds] = get_alignments(sample_1,sample_2,id)
    [score,alignment] = nwalign(nt2aa(sample_1(id).Sequence,'GeneticCode', 5),nt2aa(sample_2(id).Sequence,'GeneticCode', 5));
    x = seqinsertgaps(sample_1(id).Sequence,alignment(1,:));
    y = seqinsertgaps(sample_2(id).Sequence,alignment(3,:));
    [dn, ds, vardn, vards] = dnds(x, y, 'window', 20);
end