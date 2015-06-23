fid = fopen([ 'DbFitAbs_YYYYMMDD_augNNN_MMDDYY.m' ], 'r');
tline = fgets(fid);
tlines = '';
while ischar(tline)
    tlines = [tlines tline];
    tline = fgets(fid);
end
fclose(fid);
filesplit = strsplit(tlines,'##');
output = '';
params = readtable([ '../studydata/AUG095_032415params.csv' ],'ReadVariableNames',0,'ReadRowNames',1);


for d = 1
    p = params{:,d};
    p_format(1) = p(1);
    p_format(2) = p(10);
    p_format(3) = p(11);
    p_format(4) = p(3);
    p_format(5) = p(12);
    rN = [0 0 0 0]; %L R Lp Rp
    regionlabels = {'Left Forehead','Right Forehead','Left Parietal','Right Parietal'};
    for r = 1:4
        if length(p{12+r})
            rN(r) = 1;
        end
    end
    rSTRcell = cellstr(regionlabels(find(rN)));
    rMARKS = '';
    for r = find(rN)
        if isempty(rMARKS)
            rMARKS = ['[' p{12+r} '],'];
        else
            rMARKS = [rMARKS '[' p{12+r} ']'];
        end
    end
    p_format(6) = {rMARKS};
    rSTR = ['''' rSTRcell{1} ''''];
    for rs = 2:length(rSTRcell)
        rSTR = [rSTR ',''' rSTRcell(rs) ''''];
    end
    p_format(7) = {[rSTR{:}]};
    output = filesplit{1};
    for f = 2:length(filesplit)
        output = [output p_format{f-1} filesplit{f}];
    end
    fopen(sprintf('DbFitAbs_%s_aug%s_MMDDYY.m.m', 'w');
    w1 = fopen('test.m', 'w');
    fwrite(w1,output)
    fclose(w1)

end

