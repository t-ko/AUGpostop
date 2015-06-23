rownames = {};
mua_list = [];
muastd_list = [];
musp_list = [];
muspstd_list = [];
for rn = 1:length(grouplabels)
    if ~isempty(findstr(grouplabels{rn},'eft')) && ~isempty(findstr(grouplabels{rn},'orehead'))
        mua_list = [mua_list; nanmean(leftmua,1)];
        muastd_list = [muastd_list; nanstd(leftmua,1)];
        musp_list = [musp_list; nanmean(leftmusp,1)];
        muspstd_list = [muspstd_list; nanstd(leftmusp,1)];
        rownames = {rownames{:} grouplabels{rn}};
    elseif ~isempty(findstr(grouplabels{rn},'eft')) && (~isempty(findstr(grouplabels{rn},'ietal')) || ~isempty(findstr(grouplabels{rn},'ide')))
        mua_list = [mua_list; nanmean(leftparietalmua,1)];
        muastd_list = [muastd_list; nanstd(leftparietalmua,1)];
        musp_list = [musp_list; nanmean(leftparietalmusp,1)];
        muspstd_list = [muspstd_list; nanstd(leftparietalmusp,1)];
        rownames = {rownames{:} grouplabels{rn}};
    elseif ~isempty(findstr(grouplabels{rn},'ight')) && ~isempty(findstr(grouplabels{rn},'orehead'))
        mua_list = [mua_list; nanmean(rightmua,1)];
        muastd_list = [muastd_list; nanstd(rightmua,1)];
        musp_list = [musp_list; nanmean(rightmusp,1)];
        muspstd_list = [muspstd_list; nanstd(rightmusp,1)];
        rownames = {rownames{:} grouplabels{rn}};
    elseif ~isempty(findstr(grouplabels{rn},'ight')) && (~isempty(findstr(grouplabels{rn},'ietal')) || ~isempty(findstr(grouplabels{rn},'ide')))
        mua_list = [mua_list; nanmean(righparietaltmua,1)];
        muastd_list = [muastd_list; nanstd(rightparietalmua,1)];
        musp_list = [musp_list; nanmean(rightparietalmusp,1)];
        muspstd_list = [muspstd_list; nanstd(rightparietalmusp,1)];
        rownames = {rownames{:} grouplabels{rn}};
    end
end    

table(mua_list,muastd_list ,'RowNames',rownames,'VariableNames',{'mua', 'muaSTD'})
table(musp_list,muspstd_list ,'RowNames',rownames,'VariableNames',{'musp', 'muspSTD'})
