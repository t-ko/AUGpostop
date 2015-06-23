%%Time-stamp: "2007-05-15 15:14:13 matlabuser"

close all
clear all

studyID = 'aug068_parietal';
analysisIDs={'postop_2';'postop_6'}; %R1, R2, R3
plotID='absprops';

StO2left = NaN(1,length(analysisIDs));
StO2left_std = NaN(1,length(analysisIDs));
StO2right = NaN(1,length(analysisIDs));
StO2right_std = NaN(1,length(analysisIDs));
StO2leftparietal = NaN(1,length(analysisIDs));
StO2leftparietal_std = NaN(1,length(analysisIDs));
StO2rightparietal = NaN(1,length(analysisIDs));
StO2rightparietal_std = NaN(1,length(analysisIDs));

THCleft = NaN(1,length(analysisIDs));
THCleft_std = NaN(1,length(analysisIDs));
THCright = NaN(1,length(analysisIDs));
THCright_std = NaN(1,length(analysisIDs));
THCleftparietal = NaN(1,length(analysisIDs));
THCleftparietal_std = NaN(1,length(analysisIDs));
THCrightparietal = NaN(1,length(analysisIDs));
THCrightparietal_std = NaN(1,length(analysisIDs));

savefigures=1;
load('colors.mat');

for R = 1:length(analysisIDs)
    analysisID = analysisIDs{R};
    ext='';
    fname1=[ studyID '_' analysisID ext '_baselines.mat'];
    load(fname1); 

    fdir=['..\' studyID '\' analysisID '\'];
    files=dir([ fdir 'Data_*.txt']);
    fname=files(1).name(6:end-4);
    DateVector(R,:) = datevec([fname(1:2) '/' fname(3:4) '/' fname(5:6) ' ' fname(8:9) ':' fname(10:11) ':' fname(12:13)]);
    
    %% Process individual region trends    
    if exist('StO2_left','var') 
        StO2left(R) = nanmean(StO2_left);
        StO2left_std(R) = nanstd(StO2_left);
        THCleft(R) = nanmean(THC_left);
        THCleft_std(R) = nanstd(THC_left);
    end
    if exist('StO2_leftparietal','var') 
        StO2leftparietal(R) = nanmean(StO2_leftparietal);            
        StO2leftparietal_std(R) = nanstd(StO2_leftparietal);  
        THCleftparietal(R) = nanmean(THC_leftparietal);            
        THCleftparietal_std(R) = nanstd(THC_leftparietal); 
    end
    if exist('StO2_right','var')       
        StO2right(R) = nanmean(StO2_right);            
        StO2right_std(R) = nanstd(StO2_right); 
        THCright(R) = nanmean(THC_right);            
        THCright_std(R) = nanstd(THC_right);
    end
    if exist('StO2_rightparietal','var') 
        StO2rightparietal(R) = nanmean(StO2_rightparietal);            
        StO2rightparietal_std(R) = nanstd(StO2_rightparietal); 
        THCrightparietal(R) = nanmean(THC_rightparietal);            
        THCrightparietal_std(R) = nanstd(THC_rightparietal);  
    end
      
      
end

%% Extract post-operative time
import_AUG_surgery;
% Operative time from off X-clamp - estimate at noon if not recorded
tempdata = strsplit(studyID,'_'); 
idmatch = regexpi(StudyID,upper(tempdata{1}));
for idn = 1:length(idmatch)
    if ~isempty(idmatch{idn})
        break
    end
end
if isnan(DateofSurgery(idn))
    DateofSurgery(idn) = [DateVector(1,1:3) 0 0 0];
end
if isnan(AnesthesiaRecordTimeoffXClamp(idn))
    AnesthesiaRecordTimeoffXClamp(idn) = datenum([0 0 0 12 0 0]);
end

SurgeryVector = datevec(DateofSurgery(idn))+datevec(AnesthesiaRecordTimeoffXClamp(idn));

t_postop = [];
for p = 1:length(analysisIDs)
    t_postop(p) = datenum(DateVector(p,:)-(SurgeryVector)); %days since surgery
end
%%%%%

%%%%%%%%%%
%% PLOT TREND StO2
f=figure(1);
regionlabels = {'Left Forehead','Right Forehead','Left Parietal','Right Parietal'};
% R = 1:length(analysisIDs);
R = t_postop;
hold on;
plot(R,StO2left,'MarkerSize',20,'LineWidth',3,'Color',colors(1,:))
plot(R,StO2right,'MarkerSize',20,'LineWidth',3,'Color',colors(2,:))
plot(R,StO2leftparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(3,:))
plot(R,StO2rightparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(4,:))
legend(regionlabels);
axis tight
set(gca,'FontSize',24)
ylabel('StO_2 (%)')
xlabel('Post-Operative Day')
tmplim=get(gca,'YLim'); 
grid on
set(f,'PaperPositionMode','Auto')
maxwindows(f);
if savefigures
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_regionflow_StO2.fig'],'fig')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_regionflow_StO2.eps'],'epsc2')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_regionflow_StO2.png'],'png')
end

%% Process average bilateral trend
StO2forehead = nanmean([StO2left; StO2right]);
StO2forehead_std = sqrt(StO2left_std.^2+StO2right_std.^2);
StO2parietal = nanmean([StO2leftparietal; StO2rightparietal]);
StO2parietal_std = sqrt(StO2leftparietal_std.^2+StO2rightparietal_std.^2);

rStO2forehead = StO2forehead./StO2forehead(min(find(~isnan(StO2forehead)))).*100;
rStO2parietal = StO2parietal./StO2forehead(min(find(~isnan(StO2forehead)))).*100;

f=figure(2);
regionlabels = {'Forehead','Parietal'};
% R = 1:length(analysisIDs);
hold on;
plot(R,StO2forehead,'MarkerSize',20,'LineWidth',3,'Color',colors(1,:))
plot(R,StO2parietal,'MarkerSize',20,'LineWidth',3,'Color',colors(2,:))
legend(regionlabels);
axis tight
set(gca,'FontSize',24)
ylabel('StO_2 (%)')
xlabel('Post-Operative Day')
tmplim=get(gca,'YLim'); 
grid on
set(f,'PaperPositionMode','Auto')
maxwindows(f);
if savefigures
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_StO2.fig'],'fig')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_StO2.eps'],'epsc2')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_StO2.png'],'png')
end

f=figure(3);
regionlabels = {'Forehead','Parietal'};
% R = 1:length(analysisIDs);
hold on;
plot(R,rStO2forehead,'MarkerSize',20,'LineWidth',3,'Color',colors(1,:))
plot(R,rStO2parietal,'MarkerSize',20,'LineWidth',3,'Color',colors(2,:))
legend(regionlabels);
axis tight
set(gca,'FontSize',24)
ylabel('rStO_2 (%)')
xlabel('Post-Operative Day')
tmplim=get(gca,'YLim'); 
grid on
set(f,'PaperPositionMode','Auto')
maxwindows(f);
if savefigures
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_StO2.fig'],'fig')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_StO2.eps'],'epsc2')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_StO2.png'],'png')
end

%%%%%%%%%%
%% PLOT TREND THC
f=figure(4);
regionlabels = {'Left Forehead','Right Forehead','Left Parietal','Right Parietal'};
% R = 1:length(analysisIDs);
hold on;
plot(R,THCleft,'MarkerSize',20,'LineWidth',3,'Color',colors(1,:))
plot(R,THCright,'MarkerSize',20,'LineWidth',3,'Color',colors(2,:))
plot(R,THCleftparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(3,:))
plot(R,THCrightparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(4,:))
legend(regionlabels);
axis tight
set(gca,'FontSize',24)
ylabel('THC')
xlabel('Post-Operative Day')
tmplim=get(gca,'YLim'); 
grid on
set(f,'PaperPositionMode','Auto')
maxwindows(f);
if savefigures
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_regionflow_THC.fig'],'fig')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_regionflow_THC.eps'],'epsc2')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_regionflow_THC.png'],'png')
end

%% Process average bilateral trend
THCforehead = nanmean([THCleft; THCright]);
THCforehead_std = sqrt(THCleft_std.^2+THCright_std.^2);
THCparietal = nanmean([THCleftparietal; THCrightparietal]);
THCparietal_std = sqrt(THCleftparietal_std.^2+THCrightparietal_std.^2);

rTHCforehead = THCforehead./THCforehead(min(find(~isnan(THCforehead)))).*100;
rTHCparietal = THCparietal./THCforehead(min(find(~isnan(THCforehead)))).*100;

f=figure(5);
regionlabels = {'Forehead','Parietal'};
% R = 1:length(analysisIDs);
hold on;
plot(R,THCforehead,'MarkerSize',20,'LineWidth',3,'Color',colors(1,:))
plot(R,THCparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(2,:))
legend(regionlabels);
axis tight
set(gca,'FontSize',24)
ylabel('THC')
xlabel('Post-Operative Day')
tmplim=get(gca,'YLim'); 
grid on
set(f,'PaperPositionMode','Auto')
maxwindows(f);
if savefigures
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_THC.fig'],'fig')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_THC.eps'],'epsc2')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_THC.png'],'png')
end

f=figure(6);
regionlabels = {'Forehead','Parietal'};
% R = 1:length(analysisIDs);
hold on;
plot(R,rTHCforehead,'MarkerSize',20,'LineWidth',3,'Color',colors(1,:))
plot(R,rTHCparietal,'MarkerSize',20,'LineWidth',3,'Color',colors(2,:))
legend(regionlabels);
axis tight
set(gca,'FontSize',24)
ylabel('rTHC (%)')
xlabel('Post-Operative Day')
tmplim=get(gca,'YLim'); 
grid on
set(f,'PaperPositionMode','Auto')
maxwindows(f);
if savefigures
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_THC.fig'],'fig')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_THC.eps'],'epsc2')
    saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_globalflow_THC.png'],'png')
end

%%%%%%%%%%


tt=['save plotAbsProps_' studyID '_' plotID '.mat studyID analysisIDs t_postop StO2left StO2left_std StO2right StO2right_std StO2leftparietal StO2leftparietal_std StO2rightparietal StO2rightparietal_std StO2forehead StO2forehead_std StO2parietal StO2parietal_std rStO2forehead rStO2parietal THCleft THCleft_std THCright THCright_std THCleftparietal THCleftparietal_std THCrightparietal THCrightparietal_std THCforehead THCforehead_std THCparietal THCparietal_std rTHCforehead rTHCparietal'];
eval(tt);
