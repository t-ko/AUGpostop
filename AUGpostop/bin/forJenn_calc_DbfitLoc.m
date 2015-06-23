for rn = 1:length(regionmarks)
    rmarks = Marksflow(regionmarks{rn});
    if ~isempty(findstr(regionlabels{rn},'eft')) && ~isempty(findstr(regionlabels{rn},'orehead'))
        DbfitLoc.left = nanmean(nanmean(Dbfitavg(Marksflow(regionmarks{rn}(1)):Marksflow(regionmarks{rn}(2)))));
        DbfitLoc.left_std = nanstd(Dbfitavg(Marksflow(regionmarks{rn}(1)):Marksflow(regionmarks{rn}(2))));
    elseif ~isempty(findstr(regionlabels{rn},'eft')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
        DbfitLoc.leftparietal = nanmean(Dbfitavg(Marksflow(regionmarks{rn}(1)):Marksflow(regionmarks{rn}(2))));
        DbfitLoc.leftparietal_std = nanstd(Dbfitavg(Marksflow(regionmarks{rn}(1)):Marksflow(regionmarks{rn}(2))));
    elseif ~isempty(findstr(regionlabels{rn},'ight')) && ~isempty(findstr(regionlabels{rn},'orehead'))
        DbfitLoc.right = nanmean(Dbfitavg(Marksflow(regionmarks{rn}(1)):Marksflow(regionmarks{rn}(2))));
        DbfitLoc.right_std = nanstd(Dbfitavg(Marksflow(regionmarks{rn}(1)):Marksflow(regionmarks{rn}(2))));
    elseif ~isempty(findstr(regionlabels{rn},'ight')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
        DbfitLoc.rightparietal = nanmean(Dbfitavg(Marksflow(regionmarks{rn}(1)):Marksflow(regionmarks{rn}(2))));
        DbfitLoc.rightparietal_std = nanstd(Dbfitavg(Marksflow(regionmarks{rn}(1)):Marksflow(regionmarks{rn}(2))));
    end
end    
DbfitLoc