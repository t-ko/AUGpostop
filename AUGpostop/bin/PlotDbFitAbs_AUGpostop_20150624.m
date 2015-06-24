%%Time-stamp: "2007-05-15 15:14:13 matlabuser"

close all
clear all

ALLstudyIDs={'AUG089_120514';'AUG090_121514';'AUG091_011615';'AUG092_012615';'AUG093_020615';'AUG094_021815';'AUG095_032415';'AUG096_041315';'AUG097_041315';'AUG098_042115';'AUG099_042415';'AUG100_051415';'AUG101_052815';'AUG102_061015'};
plotID_flow='BFIabs_fitavg';
plotID_oxy='absprops';

savefigures=1;
load('colors.mat');

xrange = [0 0];
rBFIyrange = [Inf 0];
markers = {'+','o','*','x','>','^'};

ALL_BFIleft = [];
ALL_BFIright = [];
ALL_BFIleftparietal = [];
ALL_BFIrightparietal = [];
ALL_StO2left = [];
ALL_StO2right = [];
ALL_StO2leftparietal = [];
ALL_StO2rightparietal = [];
ALL_THCleft = [];
ALL_THCright = [];
ALL_THCleftparietal = [];
ALL_THCrightparietal = [];

for sID = 1:length(ALLstudyIDs)
    studyID = ALLstudyIDs{sID};
    fname1=[ 'plotDbFitAbs_' studyID '_' plotID_flow '.mat'];
    load([ fname1 ]); 
    fname2=[ 'plotAbsProps_' studyID '_' plotID_oxy '.mat'];
    load([ fname2 ]); 
    ALL_BFIleft = [ALL_BFIleft BFIleft];
    ALL_BFIright = [ALL_BFIright BFIright];
    ALL_BFIleftparietal = [ALL_BFIleftparietal BFIleftparietal];
    ALL_BFIrightparietal = [ALL_BFIrightparietal BFIrightparietal];
    ALL_StO2left = [ALL_StO2left StO2left];
    ALL_StO2right = [ALL_StO2right StO2right];
    ALL_StO2leftparietal = [ALL_StO2leftparietal StO2leftparietal];
    ALL_StO2rightparietal = [ALL_StO2rightparietal StO2rightparietal];
    ALL_THCleft = [ALL_THCleft THCleft];
    ALL_THCright = [ALL_THCright THCright];
    ALL_THCleftparietal = [ALL_THCleftparietal THCleftparietal];
    ALL_THCrightparietal = [ALL_THCrightparietal THCrightparietal];
    
    R = 1:length(analysisIDs)-1;
    if max(R) > xrange(2)
        xrange(2) = max(R);
    end
end

%% REGIONAL ANALYSIS
% StO2, compare absolute values
% nanmean(ALL_StO2rightparietal)
% ind1 = find(~isnan(ALL_StO2rightparietal));
% nanmean(ALL_StO2right(ind1))
% 
% nanmean(ALL_StO2leftparietal)
% ind1 = find(~isnan(ALL_StO2leftparietal));
% nanmean(ALL_StO2left(ind1))
% 
% % THC, compare absolute values change (
% nanmean(ALL_THCrightparietal)
% ind1 = find(~isnan(ALL_THCrightparietal));
% nanmean(ALL_THCright(ind1))
% 
% nanmean(ALL_THCleftparietal)
% ind1 = find(~isnan(ALL_THCleftparietal));
% nanmean(ALL_THCleft(ind1))
% 
% % BFI, compare to left front
% nanmean(ALL_BFIrightparietal)
% ind1 = find(~isnan(ALL_BFIrightparietal));
% nanmean(ALL_BFIright(ind1))
% 
% nanmean(ALL_BFIleftparietal)
% ind1 = find(~isnan(ALL_BFIleftparietal));
% nanmean(ALL_BFIleft(ind1))

ALL_time_oxy = {};
ALL_time_flow = {};
ALL_StO2forehead = {};
ALL_THCforehead = {};
ALL_BFIforehead = {};
ALL_Hbforehead = {};
ALL_HbO2forehead = {};
ALL_muspforehead = {};
ALL_muaforehead = {};
import_AUG_surgery;
for sID = 1:length(ALLstudyIDs)
    studyID = ALLstudyIDs{sID};
    fname1=[ 'plotDbFitAbs_' studyID '_' plotID_flow '.mat'];
    load([ fname1 ]); 
    R_flow = t_postop;
    ALL_time_flow(sID) = {t_postop};
    fname2=[ 'plotAbsProps_' studyID '_' plotID_oxy '.mat'];
    load([ fname2 ]); 
    R_oxy = t_postop;
    ALL_time_oxy(sID) = {t_postop};
    ALL_StO2forehead(sID) = {StO2forehead};
    ALL_THCforehead(sID) = {THCforehead};
    ALL_Hbforehead(sID) = {nanmean([Hbleft; Hbright])};
    ALL_HbO2forehead(sID) = {nanmean([HbO2left; HbO2right])};
    ALL_muspforehead(sID) = {nanmean([muspleft; muspright])};
    ALL_muaforehead(sID) = {nanmean([mualeft; muaright])};
    ALL_BFIforehead(sID) = {BFIforehead};
    
    R = (1:length(analysisIDs))-1;  
    
    figure(1); hold on
%     plot(R,BFIleft,'-o','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     plot(R,BFIright,'-x','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     plot(R,BFIleftparietal,'-+','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     plot(R,BFIrightparietal,'-*','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     xlim(xrange);
    
    tempdata = strsplit(studyID,'_'); 
    idmatch = regexpi(StudyID,upper(tempdata{1}));
    for idn = 1:length(idmatch)
        if ~isempty(idmatch{idn})
            break
        end
    end
    legend_diagnosis(sID) = PrimaryDiagnosis(idn);

    figure(1); hold on
    ylabel('StO_2 (%)','FontSize',25)
    errorbar(R_oxy,StO2forehead,StO2forehead_std,'-o','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     plot(R,StO2parietal,'-x','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
    %% Linear Regression
    [r,m,b] = regression(R_oxy,StO2forehead);
    slope(sID) = m;
    offset(sID) = b;
    
    figure(2); hold on
    ylabel(texlabel('Total Hemoglobin Conc. ({mu}mol/L)'),'FontSize',25)
    errorbar(R_oxy,THCforehead,THCforehead_std,'-o','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     plot(R,THCparietal,'-x','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))

    figure(3); hold on
    ylabel('Blood Flow Index (cm^2/s)','FontSize',25)
    errorbar(R_flow,BFIforehead,BFIforehead_std,'-o','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     plot(R,BFIparietal,'-x','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))

    figure(4); hold on
    [AX,h1,h2] = plotyy(R_oxy,rStO2forehead-100,R_oxy,rBFIforehead-100);
    set(h1,'LineStyle','-','LineWidth',3,'MarkerSize',10)
    set(h2,'LineStyle','--','LineWidth',3,'MarkerSize',10)    
    set(AX(1),'XLim',xrange,'FontSize',25,'YTick',-100:50:200,'YLim',[-100 200])
    set(AX(2),'XLim',xrange,'YTick',-100:50:200,'YLim',[-100 200],'FontSize',25)
    set(get(AX(1),'Ylabel'),'String','rStO_2 (%)','FontSize',25)
    set(get(AX(2),'Ylabel'),'String','rBFI (%)','FontSize',25)
    
    figure(5); hold on
    [AX,h1,h2] = plotyy(R,rBFIforehead,R,rStO2forehead);
    set(h1,'LineStyle','-','MarkerSize',10,'LineWidth',3) 
    set(h2,'LineStyle','--','MarkerSize',10,'LineWidth',3) 
end    


for fn = 1:2
    f = figure(fn);
    xlabel('Post-Operative Day','FontSize',25)
    set(gca,'FontSize',25)
    legend(legend_diagnosis,'FontSize',20)
    grid on
    set(f,'PaperPositionMode','Auto')
    maxwindows(f);
    if savefigures
        saveas(gcf,['savedfigs' filesep 'Abs_' studyID '_' plotID_oxy '_plot' num2str(fn) '.fig'],'fig')
        saveas(gcf,['savedfigs' filesep 'Abs_' studyID '_' plotID_oxy '_plot' num2str(fn) '.eps'],'epsc2')
        saveas(gcf,['savedfigs' filesep 'Abs_' studyID '_' plotID_oxy '_plot' num2str(fn) '.png'],'png')
    end
end

%%%%%%%%%%

% ALL_time_oxy = {};
% ALL_time_flow = {};
% ALL_StO2forehead = {};
% ALL_THCforehead = {};
% ALL_BFIforehead = {};
% ALL_Hbforehead = {};
% ALL_HbO2forehead = {};
% ALL_muspforehead = {};


% tt=['save plotDbFitAbs_' studyID '_' plotID '.mat names avgfit studyIDs analysisIDs sources baselinemarks markstoshow markstolabel excludemarks t0 frameperiod region fMarks BFIleft BFIleft_std BFIright BFIright_std BFIleftparietal BFIleftparietal_std BFIrightparietal BFIrightparietal_std BFIforehead BFIforehead_std BFIparietal BFIparietal_std rBFIforehead rBFIparietal'];
% eval(tt);
