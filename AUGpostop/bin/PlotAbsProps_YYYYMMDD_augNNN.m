%%Time-stamp: "2007-05-15 15:14:13 matlabuser"

close all
clear all

studyID = '##';
analysisIDs={##}; %R1, R2, R3
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

muspleft = NaN(1,length(analysisIDs));
muspleft_std = NaN(1,length(analysisIDs));
muspright = NaN(1,length(analysisIDs));
muspright_std = NaN(1,length(analysisIDs));
muspleftparietal = NaN(1,length(analysisIDs));
muspleftparietal_std = NaN(1,length(analysisIDs));
musprightparietal = NaN(1,length(analysisIDs));
musprightparietal_std = NaN(1,length(analysisIDs));

mualeft = NaN(1,length(analysisIDs));
mualeft_std = NaN(1,length(analysisIDs));
muaright = NaN(1,length(analysisIDs));
muaright_std = NaN(1,length(analysisIDs));
mualeftparietal = NaN(1,length(analysisIDs));
mualeftparietal_std = NaN(1,length(analysisIDs));
muarightparietal = NaN(1,length(analysisIDs));
muarightparietal_std = NaN(1,length(analysisIDs));

Hbleft = NaN(1,length(analysisIDs));
Hbleft_std = NaN(1,length(analysisIDs));
Hbright = NaN(1,length(analysisIDs));
Hbright_std = NaN(1,length(analysisIDs));
Hbleftparietal = NaN(1,length(analysisIDs));
Hbleftparietal_std = NaN(1,length(analysisIDs));
Hbrightparietal = NaN(1,length(analysisIDs));
Hbrightparietal_std = NaN(1,length(analysisIDs));

HbO2left = NaN(1,length(analysisIDs));
HbO2left_std = NaN(1,length(analysisIDs));
HbO2right = NaN(1,length(analysisIDs));
HbO2right_std = NaN(1,length(analysisIDs));
HbO2leftparietal = NaN(1,length(analysisIDs));
HbO2leftparietal_std = NaN(1,length(analysisIDs));
HbO2rightparietal = NaN(1,length(analysisIDs));
HbO2rightparietal_std = NaN(1,length(analysisIDs));

savefigures=1;
load('colors.mat');

for R = 1:length(analysisIDs)
    analysisID = analysisIDs{R};
    ext='';
    fname1=[ studyID '_' analysisID ext '_baselines.mat'];
    load(fname1); 

    fdir=['..' filesep '' studyID filesep analysisID filesep];
    files=dir([ fdir 'Data_*.txt']);
    fname=files(1).name(6:end-4);
    if length(fname)>13
        DateVector(R,:) = datevec([fname(5:6) '/' fname(7:8) '/' fname(1:4) ' ' fname(end-5:end-4) ':' fname(end-3:end-2) ':' fname(end-1:end)]);
    else
        DateVector(R,:) = datevec([fname(1:2) '/' fname(3:4) '/' fname(5:6) ' ' fname(end-5:end-4) ':' fname(end-3:end-2) ':' fname(end-1:end)]);
    end   
    %% Process individual region trends    
    if exist('StO2_left','var') 
        StO2left(R) = nanmean(StO2_left);
        StO2left_std(R) = nanstd(StO2_left);
        THCleft(R) = nanmean(THC_left);
        THCleft_std(R) = nanstd(THC_left);
        Hbleft(R) = nanmean(Hb_left);
        Hbleft_std(R) = nanstd(Hb_left);
        HbO2left(R) = nanmean(HbO2_left);
        HbO2left_std(R) = nanstd(HbO2_left);
        temp = nanmean(leftmusp,1);
        muspleft(R) = temp(2); %785nm
        temp = nanstd(leftmusp,1);
        muspleft_std(R) = temp(2);
        temp = nanmean(leftmua,1);
        mualeft(R) = temp(2); %785nm
        temp = nanstd(leftmua,1);
        mualeft_std(R) = temp(2);
        
    end
    if exist('StO2_leftparietal','var') 
        StO2leftparietal(R) = nanmean(StO2_leftparietal);            
        StO2leftparietal_std(R) = nanstd(StO2_leftparietal);  
        THCleftparietal(R) = nanmean(THC_leftparietal);            
        THCleftparietal_std(R) = nanstd(THC_leftparietal); 
        muspleftparietal(R) = nanmean(leftparietalmusp(:,2));
        muspleftparietal_std(R) = nanstd(leftparietalmusp(:,2));
        mualeftparietal(R) = nanmean(leftparietalmua(:,2));
        mualeftparietal_std(R) = nanstd(leftparietalmua(:,2));
    end
    if exist('StO2_right','var')       
        StO2right(R) = nanmean(StO2_right);            
        StO2right_std(R) = nanstd(StO2_right); 
        THCright(R) = nanmean(THC_right);            
        THCright_std(R) = nanstd(THC_right);
        Hbright(R) = nanmean(Hb_right);
        Hbright_std(R) = nanstd(Hb_right);
        HbO2right(R) = nanmean(HbO2_right);
        HbO2right_std(R) = nanstd(HbO2_right);
        temp = nanmean(rightmusp,1);
        muspright(R) = temp(2);
        temp = nanstd(rightmusp,1);
        muspright_std(R) = temp(2);
        temp = nanmean(rightmua,1);
        muaright(R) = temp(2);
        temp = nanstd(rightmua,1);
        muaright_std(R) = temp(2);
    end
    if exist('StO2_rightparietal','var') 
        StO2rightparietal(R) = nanmean(StO2_rightparietal);            
        StO2rightparietal_std(R) = nanstd(StO2_rightparietal); 
        THCrightparietal(R) = nanmean(THC_rightparietal);            
        THCrightparietal_std(R) = nanstd(THC_rightparietal);  
        musprightparietal(R) = nanmean(rightparietalmusp(:,2));
        musprightparietal_std(R) = nanstd(rightparietalmusp(:,2));
        muarightparietal(R) = nanmean(rightparietalmua(:,2));
        muarightparietal_std(R) = nanstd(rightparietalmua(:,2));
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
if isnan(DateofSurgery(idn)) || (sum(datevec(DateofSurgery(idn))<([2014 0 0 0 0 0]))>0)
    DateofSurgery(idn) = datenum([DateVector(1,1:3) 0 0 0]);
end
if isnan(AnesthesiaRecordTimeoffXClamp(idn))
    AnesthesiaRecordTimeoffXClamp(idn) = datenum([0 0 0 DateVector(1,4)-1 DateVector(1,5:6)]);
end

SurgeryVector = datevec(DateofSurgery(idn))+datevec(AnesthesiaRecordTimeoffXClamp(idn));

t_postop = [];
for p = 1:length(analysisIDs)
    t_postop(p) = datenum(DateVector(p,:))-datenum((SurgeryVector)); %days since surgery
end
%%%%%

%%%%%%%%%%
if savefigures
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
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_region_StO2.fig'],'fig')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_region_StO2.eps'],'epsc2')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_region_StO2.png'],'png')
    end
end

%% Process average bilateral trend
StO2forehead = nanmean([StO2left; StO2right]);
StO2forehead_std = sqrt(StO2left_std.^2+StO2right_std.^2);
StO2parietal = nanmean([StO2leftparietal; StO2rightparietal]);
StO2parietal_std = sqrt(StO2leftparietal_std.^2+StO2rightparietal_std.^2);

rStO2forehead = StO2forehead./StO2forehead(min(find(~isnan(StO2forehead)))).*100;
rStO2parietal = StO2parietal./StO2forehead(min(find(~isnan(StO2forehead)))).*100;

if savefigures
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
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_StO2.fig'],'fig')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_StO2.eps'],'epsc2')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_StO2.png'],'png')
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
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_StO2.fig'],'fig')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_StO2.eps'],'epsc2')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_StO2.png'],'png')
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
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_region_THC.fig'],'fig')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_region_THC.eps'],'epsc2')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_region_THC.png'],'png')
    end
end

%% Process average bilateral trend
THCforehead = nanmean([THCleft; THCright]);
THCforehead_std = sqrt(THCleft_std.^2+THCright_std.^2);
THCparietal = nanmean([THCleftparietal; THCrightparietal]);
THCparietal_std = sqrt(THCleftparietal_std.^2+THCrightparietal_std.^2);

rTHCforehead = THCforehead./THCforehead(min(find(~isnan(THCforehead)))).*100;
rTHCparietal = THCparietal./THCforehead(min(find(~isnan(THCforehead)))).*100;

if savefigures
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
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_THC.fig'],'fig')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_THC.eps'],'epsc2')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_THC.png'],'png')
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
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_THC.fig'],'fig')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_THC.eps'],'epsc2')
        saveas(gcf,['../' studyID '/savedfigs/AbsProps_' studyID '_' plotID '_global_THC.png'],'png')
    end
end

%%%%%%%%%%


tt=['save plotAbsProps_' studyID '_' plotID '.mat studyID analysisIDs t_postop StO2left StO2left_std StO2right StO2right_std StO2leftparietal StO2leftparietal_std StO2rightparietal StO2rightparietal_std StO2forehead StO2forehead_std StO2parietal StO2parietal_std rStO2forehead rStO2parietal THCleft THCleft_std THCright THCright_std THCleftparietal THCleftparietal_std THCrightparietal THCrightparietal_std THCforehead THCforehead_std THCparietal THCparietal_std rTHCforehead rTHCparietal muspleft muspleft_std muspright muspright_std muspleftparietal muspleftparietal_std musprightparietal musprightparietal_std mualeft mualeft_std muaright muaright_std mualeftparietal mualeftparietal_std muarightparietal muarightparietal_std Hbleft Hbleft_std Hbright Hbright_std Hbleftparietal Hbleftparietal_std Hbrightparietal Hbrightparietal_std HbO2left HbO2left_std HbO2right HbO2right_std HbO2leftparietal HbO2leftparietal_std HbO2rightparietal HbO2rightparietal_std '];
eval(tt);
