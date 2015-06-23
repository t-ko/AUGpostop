function cmd = AUG_readparams(filename,func,eval_functions)
    if (nargin < 2)
        if ~isempty(strfind(filename,'DCSonly'))
            func = {'DbFitAbs','Checkfits','PlotDbFitAbs'};
        else
            func = {};
        end
    end
    temp = char(reshape(func, 1,length(func)))';
    func_str = temp(:)';
    if (nargin < 3)
        eval_functions = false;
    end
    delimiter = '##';
    params = readtable([ '../studydata/' filename ],'ReadVariableNames',0,'ReadRowNames',1);
    regionlabels = {'Left Forehead','Right Forehead','Left Parietal','Right Parietal'};
    cmd = '';
    for d = 1:size(params,2)
        p = params{:,d};
        if iscell(p(1))
            if isempty(func) || ~isempty(strfind(func_str,'OptProps'))
                %% OptProps
                filesplit = importtemplate('OptProps_YYYYMMDD_augNNN_MMDDYY.m',delimiter);
                p_format(1) = p(1);
                p_format(2) = p(2);
                p_format(3) = p(3);
                rN = [0 0 0 0]; %L R Lp Rp
                % Determine which regions have data marks
                for r = 1:4
                    if length(p{5+r})
                        rN(r) = 1;
                    end
                end
                % Populate groups string with indices
                rMARKS = ['[' p{4} '];[' p{5} ']'];
                for r = find(rN)
                    rMARKS = [rMARKS ';[' p{5+r} ']'];
                end
                p_format(4) = {rMARKS};
                % Populate grouplabel string
                rSTRcell = cellstr(regionlabels(find(rN)));
                rSTR = ['''' rSTRcell{1} ''''];
                for rs = 2:length(rSTRcell)
                    rSTR = [rSTR ';''' rSTRcell(rs) ''''];
                end
                p_format(5) = {['''Calibration block'';''Check block'';' rSTR{:}]};
                % Concatenate output strings into full file
                output = filesplit{1};
                for f = 2:length(filesplit)
                    output = [output p_format{f-1} filesplit{f}];
                end
                % Write output to formatted filename
                if ~isempty(strmatch(p{3},'absolutes'))
                    w1 = fopen(sprintf('OptProps_%s_%s_1.m',p{2},lower(p{1})), 'w');
                else
                    w1 = fopen(sprintf('OptProps_%s_%s.m',p{2},lower(p{1})), 'w');
                end
                
                fwrite(w1,output);
                fclose(w1);
                clear p_format output filesplit
            end

            if isempty(func) || ~isempty(strfind(func_str,'absprops'))
                %% absprops
                filesplit = importtemplate('absprops_YYYYMMDD_augNNN_MMDDYY.m',delimiter);
                p_format(1) = p(1);
                p_format(2) = p(3);
                % Concatenate output strings into full file
                output = filesplit{1};
                for f = 2:length(filesplit)
                    output = [output p_format{f-1} filesplit{f}];
                end
                % Write output to formatted filename
                if ~isempty(strmatch(p{3},'absolutes'))
                    w1 = fopen(sprintf('absprops_%s_%s_1.m',p{2},lower(p{1})), 'w');
                else
                    w1 = fopen(sprintf('absprops_%s_%s.m',p{2},lower(p{1})), 'w');
                end
                fwrite(w1,output);
                fclose(w1);
                clear p_format output filesplit
            end

            if isempty(func) || (~isempty(strfind(func_str,'DbFitAbs')) && isempty(strfind(func_str,'PlotDbFitAbs'))) || (length(strfind(func_str,'DbFitAbs'))>1)
                %% DbFitAbs
                filesplit = importtemplate('DbFitAbs_YYYYMMDD_augNNN_MMDDYY.m',delimiter);
                p_format(1) = p(1);
                p_format(2) = p(10);
                p_format(3) = p(11);
                p_format(4) = p(3);
                p_format(5) = p(12);
                rN = [0 0 0 0]; %L R Lp Rp
                % Determine which regions have data marks
                for r = 1:4
                    if length(p{12+r})
                        rN(r) = 1;
                    end
                end
                % Format region marks string
                rMARKS = '';
                for r = find(rN)
                    if isempty(rMARKS)
                        rMARKS = ['[' p{12+r} ']'];
                    else
                        rMARKS = [rMARKS ',[' p{12+r} ']'];
                    end
                end
                p_format(6) = {rMARKS};
                % Populate string with present regionlabels
                rSTRcell = cellstr(regionlabels(find(rN)));
                rSTR = {['''' rSTRcell{1} '''']};
                for rs = 2:length(rSTRcell)
                    rSTR = [rSTR ',''' rSTRcell(rs) ''''];
                end
                p_format(7) = {[rSTR{:}]};
                % Concatenate output strings into full file
                output = filesplit{1};
                for f = 2:length(filesplit)
                    output = [output p_format{f-1} filesplit{f}];
                end
                % Write output to formatted filename
                if ~isempty(strmatch(p{10},'BFI_1'))
                    w1 = fopen(sprintf('DbFitAbs_%s_%s_1.m',p{2},lower(p{1})), 'w');
                else
                    w1 = fopen(sprintf('DbFitAbs_%s_%s.m',p{2},lower(p{1})), 'w');
                end
                fwrite(w1,output);
                fclose(w1);
                clear p_format output filesplit
        %         run(sprintf('DbFitAbs_%s_%s.m',p{2},lower(p{1})));
            end

            if isempty(func) || ~isempty(strfind(func_str,'Checkfits'))
                %% Checkfits
                filesplit = importtemplate('Checkfits_YYYYMMDD_augNNN_MMDDYY.m',delimiter);
                p_format(1) = p(1);
                p_format(2) = p(11);
                % Concatenate output strings into full file
                output = filesplit{1};
                for f = 2:length(filesplit)
                    output = [output p_format{f-1} filesplit{f}];
                end
                % Write output to formatted filename
                if ~isempty(strmatch(p{10},'BFI_1'))
                    w1 = fopen(sprintf('Checkfits_%s_%s_1.m',p{2},lower(p{1})), 'w');
                else
                    w1 = fopen(sprintf('Checkfits_%s_%s.m',p{2},lower(p{1})), 'w');
                end
                fwrite(w1,output);
                fclose(w1);
                clear p_format output filesplit 
            end
    %         run(sprintf('Checkfits_%s_%s.m',p{2},lower(p{1})));
            cmd = [cmd sprintf('run(''OptProps_%s_%s.m'');',p{2},lower(p{1})) sprintf('run(''absprops_%s_%s.m'');',p{2},lower(p{1}))];
        else
            params = params(:,1:(d-1));
        end
    end
    if isempty(func) || ~isempty(strfind(func_str,'PlotAbsProps'))
        %% PlotAbsProps
        filesplit = importtemplate('PlotAbsProps_YYYYMMDD_augNNN.m',delimiter);
        p_format(1) = {char(params{1,1})};    
        rSTRcell = cellstr(params{3,:});
        rSTR = ['''' rSTRcell{1} ''''];
        for rs = 2:length(rSTRcell)
            rSTR = [rSTR ';''' rSTRcell(rs) ''''];
        end
        p_format(2) = {[rSTR{:}]};
        % Concatenate output strings into full file
        output = filesplit{1};
        for f = 2:length(filesplit)
            output = [output p_format{f-1} filesplit{f}];
        end
        % Write output to formatted filename
        ID = lower(params{1,1});
        ID = ID{1};
        patdate = params{2,1};
        w1 = fopen(sprintf('PlotAbsProps_%s_%s.m',patdate{1},ID(1:6)), 'w');
        fwrite(w1,output);
        fclose(w1);
        clear p_format output filesplit 
    end
    if isempty(func) || ~isempty(strfind(func_str,'PlotDbFitAbs'))
        %% PlotDbFitAbs subject summary
        filesplit = importtemplate('PlotDbFitAbs_YYYYMMDD_augNNN.m',delimiter);
        %studyID
        p_format(1) = {char(params{1,1})};
        %BFI_prefixes
        rSTRcell = cellstr(params{10,:});
        rSTR = ['''' rSTRcell{1} ''''];
        for rs = 2:length(rSTRcell)
            rSTR = [rSTR ';''' rSTRcell(rs) ''''];
        end
        p_format(2) = {[rSTR{:}]};
        %BFI_analysisIDs
        rSTRcell = cellstr(params{11,:});
        rSTR = ['''' rSTRcell{1} ''''];
        for rs = 2:length(rSTRcell)
            rSTR = [rSTR ';''' rSTRcell(rs) ''''];
        end
        p_format(3) = {[rSTR{:}]};
        % Concatenate output strings into full file
        output = filesplit{1};
        for f = 2:length(filesplit)
            output = [output p_format{f-1} filesplit{f}];
        end
        % Write output to formatted filename
        ID = lower(params{1,1});
        ID = ID{1};
        patdate = params{2,1};
        w1 = fopen(sprintf('PlotDbFitAbs_%s_%s.m',patdate{1},ID(1:6)), 'w');
        fwrite(w1,output);
        fclose(w1);
        clear p_format output filesplit 
        
    end
    
       
    if eval_functions
        eval(cmd);
    end
end