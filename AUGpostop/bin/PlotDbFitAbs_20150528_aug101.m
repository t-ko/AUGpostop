%%Time-stamp: "2007-05-15 15:14:13 matlabuser"

close all
clear all
import_AUG_surgery;

studyID = 'AUG101_052815';
prefixes={'BFI_1';'BFI_6';'BFI_052915';'BFI_053015';'BFI_053115';'BFI_060115';'BFI_060215';'BFI_060315'};
analysisIDs={'BFI_1_abs_fitavg';'BFI_052815abs_fitavg';'BFI_052915abs_fitavg';'BFI_053015abs_fitavg';'BFI_053115abs_fitavg';'BFI_060115abs_fitavg';'BFI_060215abs_fitavg';'BFI_060315abs_fitavg'}; %R1, R2, R3
avgfit=zeros(1,length(prefixes)); %0 if averaged, channel number if not
plotID='BFIabs_fitavg';
names = prefixes;
sources=zeros(1,length(prefixes)); % source selection, 1 or 2, 0 if sole source;

baselinemarks=[1 2]; % Must refer to last studyID/analysisID

autoloadmarks=1; % Autoload marks below(TK 11-11-2014)

BFIleft = NaN(1,length(prefixes));
BFIleft_std = NaN(1,length(prefixes));
BFIright = NaN(1,length(prefixes));
BFIright_std = NaN(1,length(prefixes));
BFIleftparietal = NaN(1,length(prefixes));
BFIleftparietal_std = NaN(1,length(prefixes));
BFIrightparietal = NaN(1,length(prefixes));
BFIrightparietal_std = NaN(1,length(prefixes));


excludemarks = [];
savefigures=1;
load('colors.mat');

%% Extract post-operative time
% Extract file acquisition time
t_postop = [];
for p = 1:length(prefixes)
    prefix = prefixes{p};
    fdir1 = ['..' filesep studyID '' filesep];
    fname1 = [ prefix '_' ];
    fname=[fdir1 fname1 'flow_' sprintf('%01d',0) '.dat'];
    fid = fopen([ fname ], 'r');
    [tmpdata, count]=fscanf(fid,'%c %s',7);
    fclose(fid);
    DateVector(p,:) = datevec(tmpdata(2:end));
end

% Operative time from off X-clamp - estimate at noon if not recorded
tempdata = strsplit(studyID,'_'); 
idmatch = regexpi(StudyID,upper(tempdata{1}));
for idn = 1:length(idmatch)
    if ~isempty(idmatch{idn})
        break
    end
end
if isnan(DateofSurgery(idn)) || (sum(datevec(DateofSurgery(idn))<([2014 0 0 0 0 0]))>0)
    DateofSurgery(idn) = datenum([DateVector(1,1:3) 0 0 0]);
end
if isnan(AnesthesiaRecordTimeoffXClamp(idn))
    AnesthesiaRecordTimeoffXClamp(idn) = datenum([0 0 0 DateVector(1,4)-1 DateVector(1,5:6)]);
end

SurgeryVector = datevec(DateofSurgery(idn))+datevec(AnesthesiaRecordTimeoffXClamp(idn));

t_postop = [];
for p = 1:length(prefixes)
    t_postop(p) = datenum(DateVector(p,:))-datenum((SurgeryVector)); %days since surgery
end
%%%%%

t0 = 0;

for R = 1:length(analysisIDs)
    analysisID = analysisIDs{R};
    s = sources(R);
    plotCH=avgfit(R);
    if ~plotCH
        %% AVERAGE
        fname1=[ studyID '_S' num2str(max([s 1]))];
        if ~isempty(analysisID)
            fname1 = [fname1 '_' analysisID];
        end
        load([ fname1 '_flow_output_fitavg.mat']); 
        flowdata = Dbfitavg;
    else
        %% INDIVIDUAL PLOTS
        fname1=[ studyID '_S' num2str(max([s 1]))];
        if ~isempty(analysisID)
            fname1 = [fname1 '_' analysisID];
        end
        load([ fname1 '_flow_output_fitindiv.mat']);
        flowdata = Dbfit;
    end
    
    %% Process individual region trends    
    for rn = 1:length(regionmarks)
        rmarks = Marksflow(regionmarks{rn});
        if ~isempty(findstr(regionlabels{rn},'eft')) && ~isempty(findstr(regionlabels{rn},'orehead'))
            BFIleft(R) = nanmean(flowdata(rmarks(1):rmarks(2)));
            BFIleft_std(R) = nanstd(flowdata(rmarks(1):rmarks(2)));
        elseif ~isempty(findstr(regionlabels{rn},'eft')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
            BFIleftparietal(R) = nanmean(flowdata(rmarks(1):rmarks(2)));            
            BFIleftparietal_std(R) = nanstd(flowdata(rmarks(1):rmarks(2))); 
        elseif ~isempty(findstr(regionlabels{rn},'ight')) && ~isempty(findstr(regionlabels{rn},'orehead'))
            BFIright(R) = nanmean(flowdata(rmarks(1):rmarks(2)));            
            BFIright_std(R) = nanstd(flowdata(rmarks(1):rmarks(2)));            
        elseif ~isempty(findstr(regionlabels{rn},'ight')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
            BFIrightparietal(R) = nanmean(flowdata(rmarks(1):rmarks(2)));            
            BFIrightparietal_std(R) = nanstd(flowdata(rmarks(1):rmarks(2)));            
        end
    end
    
    % Laser 1 is acquired before Laser 2
    % Adjust time-stamp accordingly

end
if ~autoloadmarks
    markstoshow=[1 2 3 4];
    markstolabel={'Left','','Right',''};
else
    markstoshow = [regionmarks{:}];
    spacelist = {};
    for i = 1:length(regionlabels)
        spacelist = [spacelist {''}];
    end
    regionlist = [regionlabels; spacelist];
    markstolabel = regionlist(:)';
end
region = [];
frameperiod = (timeaxis_flow(2)-timeaxis_flow(1)).*60; 
fMarks=(timeaxis_flow(Marksflow(markstoshow))-timeaxis_flow(1)).*60; %seconds

%%%%%%%%%%
%% PLOT TREND
R = t_postop;

f=figure(3);
regionlabels = {'Left Forehead','Right Forehead','Left Parietal','Right Parietal'};
% R = 1:length(analysisIDs);
hold on;
plot(R,BFIleft,'MarkerSize',20,'LineWidth',3,'Color',colors(1,:))
plot(R,BFIright,'MarkerSize',20,'LineWidth',3,'Color',colors(2,:))
plot(R,BFIleftparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(3,:))
plot(R,BFIrightparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(4,:))
legend(regionlabels);
axis tight
set(gca,'FontSize',24)
ylabel('BFI')
xlabel('Post-Operative Day')
tmplim=get(gca,'YLim'); 
grid on
set(f,'PaperPositionMode','Auto')
maxwindows(f);
if savefigures
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitAbs_' studyID '_' plotID '_regionflow_BFI.fig'],'fig')
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitAbs_' studyID '_' plotID '_regionflow_BFI.eps'],'epsc2')
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitAbs_' studyID '_' plotID '_regionflow_BFI.png'],'png')
end

% Process average bilateral trend
BFIforehead = nanmean([BFIleft; BFIright]);
BFIforehead_std = sqrt(BFIleft_std.^2+BFIright_std.^2);
BFIparietal = nanmean([BFIleftparietal; BFIrightparietal]);
BFIparietal_std = sqrt(BFIleftparietal_std.^2+BFIrightparietal_std.^2);

rBFIforehead = BFIforehead./BFIforehead(min(find(~isnan(BFIforehead)))).*100;
rBFIparietal = BFIparietal./BFIforehead(min(find(~isnan(BFIforehead)))).*100;

f=figure(4);
regionlabels = {'Forehead','Parietal'};
% R = 1:length(analysisIDs);
hold on;
plot(R,rBFIforehead,'MarkerSize',20,'LineWidth',3,'Color',colors(1,:))
plot(R,rBFIparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(2,:))
legend(regionlabels);
axis tight
set(gca,'FontSize',24)
ylabel('rBFI (%)')
xlabel('Post-Operative Day')
tmplim=get(gca,'YLim'); 
grid on
set(f,'PaperPositionMode','Auto')
maxwindows(f);
if savefigures
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitAbs_' studyID '_' plotID '_globalflow_BFI.fig'],'fig')
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitAbs_' studyID '_' plotID '_globalflow_BFI.eps'],'epsc2')
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitAbs_' studyID '_' plotID '_globalflow_BFI.png'],'png')
end

f=figure(5);
regionlabels = {'Forehead','Parietal'};
% R = 1:length(analysisIDs);
hold on;
plot(R,BFIforehead,'MarkerSize',20,'LineWidth',3,'Color',colors(1,:))
plot(R,BFIparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(2,:))
legend(regionlabels);
axis tight
set(gca,'FontSize',24)
ylabel('rBFI')
xlabel('Post-Operative Day')
tmplim=get(gca,'YLim'); 
grid on
set(f,'PaperPositionMode','Auto')
maxwindows(f);
if savefigures
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitAbs_' studyID '_' plotID '_globalflow_BFI.fig'],'fig')
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitAbs_' studyID '_' plotID '_globalflow_BFI.eps'],'epsc2')
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitAbs_' studyID '_' plotID '_globalflow_BFI.png'],'png')
end

%%%%%%%%%%


tt=['save plotDbFitAbs_' studyID '_' plotID '.mat names avgfit studyID analysisIDs sources baselinemarks markstoshow markstolabel excludemarks t0 frameperiod region fMarks BFIleft BFIleft_std BFIright BFIright_std BFIleftparietal BFIleftparietal_std BFIrightparietal BFIrightparietal_std BFIforehead BFIforehead_std BFIparietal BFIparietal_std rBFIforehead rBFIparietal t_postop'];
eval(tt);
