function [selection, optionsModel] = optimguiImportOptions(unused)
%OPTIMGUIIMPORTOPTIONS GUI helper to import options structure. If nargin is zero
%   then a list dialog box is displayed showing 'default options' and other valid 
%   options structure that may exist in the workspace. If nargin is one then no 
%   list dialog box is displayed; default options is returned. This is a case when 
%   user selects 'Reset Optimtool' from File menu and we do not want the list dialog 
%    to popup. 
    
%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:47:44 $

% The value of 'unused' is not used. It is merely used to check nargin

if nargin < 1
    showDialog = true;
else
    showDialog = false;
end
optionsFieldnames = fieldnames(optimset);
% We check for five matching fields (to make sure that MATLAB optimization
% solvers get through this test)
minNumberOfOptionsToCheck = 5; % Display, MaxFunEvals, MaxIter, TolX, and TolFun.

selection = '';
optionsModel = '';
if showDialog
    whoslist =evalin('base','whos');
    names = {'default options'};
    for i = 1:length(whoslist)
        if strcmp(whoslist(i).class, 'struct') && strcmp(num2str(whoslist(i).size), '1  1')
            s = evalin('base', whoslist(i).name);
            if validOptions(s,optionsFieldnames,minNumberOfOptionsToCheck)
                names{end + 1 } = whoslist(i).name;
            end
        end
    end

    [selectionIndex, Answer] = listdlg('ListString', names, 'SelectionMode', 'Single', ...
        'ListSize', [250 200], 'Name', 'Import Optimization Options', ...
        'PromptString', 'Select an options structure to import:', ...
        'OKString', 'Import');
else
    Answer = 1;
    selectionIndex = 1;
end
% Answer == 1 means that user has pressed 'Import options' button. 
if Answer == 1
    if selectionIndex == 1  %default
        selection = 'default';
        options = optimset('Display', 'off'); % This is the default for the GUI
    else
        selection = names{selectionIndex};
        options = evalin('base', selection);
    end
    % Create Java hash table which is passed back to the GUI
    optionsModel = createHashTable(options,fieldnames(optimset),selection);
    % Also save the options structure in MATLAB workspace (appdata)
    setappdata(0,'optimTool_Options_Data',options);
end
% Reset Java hashtable for options change
resetOptimtoolHashTable('optimTool_Options_HashTable');



    
    
