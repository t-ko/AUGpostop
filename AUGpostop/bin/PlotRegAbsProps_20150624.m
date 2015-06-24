%%Time-stamp: "2007-05-15 15:14:13 matlabuser"

close all
clear all

ALLstudyIDs={'aug068_parietal';'poa038_parietal';'AUG092_012615';'AUG093_020615';'AUG094_021815';'AUG095_032415';'AUG096_041315';'AUG097_041315';'AUG098_042115';'AUG099_042415';'AUG100_051415';'AUG101_052815';'AUG102_061015'};
plotID_flow='BFIabs_fitavg';
plotID_oxy='absprops';
regionlabels = {'Right_Parietal','Right_Forehead','Left_Forehead','Left_Parietal'};
regionlabels_std = {'Right_Parietal_STD','Right_Forehead_STD','Left_Forehead_STD','Left_Parietal_STD'};

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
ALL_BFIleft_std = [];
ALL_BFIright_std = [];
ALL_BFIleftparietal_std = [];
ALL_BFIrightparietal_std = [];
ALL_StO2left_std = [];
ALL_StO2right_std = [];
ALL_StO2leftparietal_std = [];
ALL_StO2rightparietal_std = [];
ALL_THCleft_std = [];
ALL_THCright_std = [];
ALL_THCleftparietal_std = [];
ALL_THCrightparietal_std = [];

for sID = 1:length(ALLstudyIDs)
    studyID = ALLstudyIDs{sID};
    fname1=[ 'plotDbFitAbs_' studyID '_' plotID_flow '.mat'];
    load([ fname1 ]); 
    R_flow = t_postop;
    fname2=[ 'plotAbsProps_' studyID '_' plotID_oxy '.mat'];
    load([ fname2 ]); 
    R_oxy = t_postop;
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
    
    ALL_BFIleft_std = [ALL_BFIleft_std BFIleft_std];
    ALL_BFIright_std = [ALL_BFIright_std BFIright_std];
    ALL_BFIleftparietal_std = [ALL_BFIleftparietal_std BFIleftparietal_std];
    ALL_BFIrightparietal_std = [ALL_BFIrightparietal_std BFIrightparietal_std];
    ALL_StO2left_std = [ALL_StO2left_std StO2left_std];
    ALL_StO2right_std = [ALL_StO2right_std StO2right_std];
    ALL_StO2leftparietal_std = [ALL_StO2leftparietal_std StO2leftparietal_std];
    ALL_StO2rightparietal_std = [ALL_StO2rightparietal_std StO2rightparietal_std];
    ALL_THCleft_std = [ALL_THCleft_std THCleft_std];
    ALL_THCright_std = [ALL_THCright_std THCright_std];
    ALL_THCleftparietal_std = [ALL_THCleftparietal_std THCleftparietal_std];
    ALL_THCrightparietal_std = [ALL_THCrightparietal_std THCrightparietal_std];
    
    
    R = 1:length(analysisIDs)-1;
    if max(R) > xrange(2)
        xrange(2) = max(R);
    end
end
ind_StO2 = find(((~isnan(ALL_StO2left))+(~isnan(ALL_StO2right))+(~isnan(ALL_StO2rightparietal))+(~isnan(ALL_StO2leftparietal)))>2);
ind_THC = find(((~isnan(ALL_THCleft))+(~isnan(ALL_THCright))+(~isnan(ALL_THCrightparietal))+(~isnan(ALL_THCleftparietal)))>2);
ind_BFI = find(((~isnan(ALL_BFIleft))+(~isnan(ALL_BFIright))+(~isnan(ALL_BFIrightparietal))+(~isnan(ALL_BFIleftparietal)))>2);

regional_StO2=array2table([ALL_StO2rightparietal(ind_StO2); ...
ALL_StO2right(ind_StO2); ...
ALL_StO2left(ind_StO2); ...
ALL_StO2leftparietal(ind_StO2)]','VariableNames',regionlabels);
writetable(regional_StO2);

regional_StO2_std=array2table([ALL_StO2rightparietal_std(ind_StO2); ...
ALL_StO2right_std(ind_StO2); ...
ALL_StO2left_std(ind_StO2); ...
ALL_StO2leftparietal_std(ind_StO2)]','VariableNames',regionlabels_std);
writetable(regional_StO2_std);

regional_dStO2=table2array(regional_StO2); %difference from right frontal
regional_dStO2_std = table2array(regional_StO2_std);
for col = 1:size(regional_dStO2,2)
    regional_dStO2(:,col) = regional_dStO2(:,col)-(ALL_StO2right(ind_StO2)');
    regional_dStO2_std(:,col) = sqrt(regional_dStO2_std(:,col).^2 + (ALL_StO2right_std(ind_StO2)'.^2));
end
reg_avg_dStO2 = nanmean(regional_dStO2,1);
reg_avg_dStO2_std = sqrt(nansum(regional_dStO2_std.^2));
f=figure(1);
subplot(1,3,1);
bp = boxplot(regional_dStO2,'labels',{'Rp','R','L','Lp'});
set(bp(:,:),'linewidth',4);
ylabel('dStO_2 (%)','FontSize',25)
set(gca,'FontSize',25)
set(findobj(gca,'Type','text'),'FontSize',25)

regional_THC=array2table([ALL_THCrightparietal(ind_THC); ...
ALL_THCright(ind_THC); ...
ALL_THCleft(ind_THC); ...
ALL_THCleftparietal(ind_THC)]','VariableNames',regionlabels);
writetable(regional_THC);

regional_THC_std=array2table([ALL_THCrightparietal_std(ind_THC); ...
ALL_THCright_std(ind_THC); ...
ALL_THCleft_std(ind_THC); ...
ALL_THCleftparietal_std(ind_THC)]','VariableNames',regionlabels_std);
writetable(regional_THC_std);

regional_dTHC=table2array(regional_THC); % difference from right frontal
regional_dTHC_std = table2array(regional_THC_std);
for col = 1:size(regional_dTHC,2)
    regional_dTHC(:,col) = regional_dTHC(:,col)-(ALL_THCright(ind_THC)');
    regional_dTHC_std(:,col) = sqrt(regional_dTHC_std(:,col).^2 + (ALL_THCright_std(ind_THC)'.^2));
end
reg_avg_dTHC = nanmean(regional_dTHC,1);
reg_avg_dTHC_std = sqrt(nansum(regional_dTHC_std.^2));
subplot(1,3,2);
bp = boxplot(regional_dTHC,'labels',{'Rp','R','L','Lp'});
ylabel(texlabel('dTHC ({mu}mol/L)'),'FontSize',25)
set(bp(:,:),'linewidth',4);
set(gca,'FontSize',25)
set(findobj(gca,'Type','text'),'FontSize',25)

regional_BFI=array2table([ALL_BFIrightparietal(ind_BFI); ...
ALL_BFIright(ind_BFI); ...
ALL_BFIleft(ind_BFI); ...
ALL_BFIleftparietal(ind_BFI)]','VariableNames',regionlabels);
writetable(regional_BFI);

regional_BFI_std=array2table([ALL_BFIrightparietal_std(ind_BFI); ...
ALL_BFIright_std(ind_BFI); ...
ALL_BFIleft_std(ind_BFI); ...
ALL_BFIleftparietal_std(ind_BFI)]','VariableNames',regionlabels_std);
writetable(regional_BFI_std);

regional_rBFI=table2array(regional_BFI); % Percentage change from right frontal
regional_rBFI_std = zeros(size(regional_rBFI));
for col = 1:size(regional_rBFI,2)
    regional_rBFI(:,col) = (regional_rBFI(:,col)./(ALL_BFIright(ind_BFI)')-1)*100;
    regional_rBFI_std(:,col) = regional_rBFI(:,col).*sqrt((regional_rBFI_std(:,col)./regional_rBFI(:,col)).^2+((ALL_BFIright_std(ind_BFI)')./(ALL_BFIright(ind_BFI)')).^2);
end
reg_avg_rBFI = nanmean(regional_rBFI,1);
reg_avg_rBFI_std = sqrt(nansum(regional_rBFI_std.^2));
subplot(1,3,3);
bp = boxplot(regional_rBFI,'labels',{'Rp','R','L','Lp'});
set(bp(:,:),'linewidth',4);
set(gca,'FontSize',25)
set(findobj(gca,'Type','text'),'FontSize',25)
ylabel('rBFI (%)','FontSize',25)


%% REGIONAL ANALYSIS
% dStO2, compare absolute values
% nanmean(ALL_StO2rightparietal)
% ind1 = find(~isnan(ALL_StO2rightparietal));
% nanmean(ALL_StO2right(ind1))
% 
% nanmean(ALL_StO2leftparietal)
% ind1 = find(~isnan(ALL_StO2leftparietal));
% nanmean(ALL_StO2left(ind1))
% 
% dTHC, compare absolute values change 
% nanmean(ALL_THCrightparietal)
% ind1 = find(~isnan(ALL_THCrightparietal));
% nanmean(ALL_THCright(ind1))
% 
% nanmean(ALL_THCleftparietal)
% ind1 = find(~isnan(ALL_THCleftparietal));
% nanmean(ALL_THCleft(ind1))
% 
% rBFI, compare to left front
% nanmean(ALL_BFIrightparietal)
% ind1 = find(~isnan(ALL_BFIrightparietal));
% nanmean(ALL_BFIright(ind1))
% 
% nanmean(ALL_BFIleftparietal)
% ind1 = find(~isnan(ALL_BFIleftparietal));
% nanmean(ALL_BFIleft(ind1))


% for sID = 1:length(ALLstudyIDs)
%     studyID = ALLstudyIDs{sID};
%     fname1=[ 'plotDbFitAbs_' studyID '_' plotID_flow '.mat'];
%     load([ fname1 ]); 
%     fname2=[ 'plotAbsProps_' studyID '_' plotID_oxy '.mat'];
%     load([ fname2 ]); 
%     
%     R = (1:length(analysisIDs))-1;  
%     
%     figure(1); hold on
%     plot(R,BFIleft,'-o','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     plot(R,BFIright,'-x','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     plot(R,BFIleftparietal,'-+','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     plot(R,BFIrightparietal,'-*','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     xlim(xrange);
%     

%     figure(1); hold on
%     ylabel('StO_2 (%)','FontSize',25)
%     plot(R,StO2forehead,'-o','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     %plot(R,StO2parietal,'-x','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
% 
%     figure(2); hold on
%     ylabel(texlabel('Total Hemoglobin Conc. ({mu}mol/L)'),'FontSize',25)
%     plot(R,THCforehead,'-o','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     %plot(R,THCparietal,'-x','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
% 
%     figure(3); hold on
%     ylabel('Blood Flow Index (cm^2/s)','FontSize',25)
%     plot(R,BFIforehead,'-o','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
%     %plot(R,BFIparietal,'-x','MarkerSize',10,'LineWidth',3,'Color',colors(sID,:))
% 
%     figure(4); hold on
%     [AX,h1,h2] = plotyy(R,rStO2forehead-100,R,rBFIforehead-100);
%     set(h1,'LineStyle','-','LineWidth',3,'MarkerSize',10)
%     set(h2,'LineStyle','--','LineWidth',3,'MarkerSize',10)    
%     set(AX(1),'XLim',xrange,'FontSize',25,'YTick',-100:50:200,'YLim',[-100 200])
%     set(AX(2),'XLim',xrange,'YTick',-100:50:200,'YLim',[-100 200],'FontSize',25)
%     set(get(AX(1),'Ylabel'),'String','rStO_2 (%)','FontSize',25)
%     set(get(AX(2),'Ylabel'),'String','rBFI (%)','FontSize',25)
%     %plotyy(R,rBFIparietal,R,rStO2forehead,'MarkerSize',10,'LineWidth',3,'Color',colors(sID,:)) 
% end    
% 
% 
% for fn = 1:4
%     f = figure(fn);
%     xlabel('Post-Operative Day','FontSize',25)
%     set(gca,'FontSize',25)
%     grid on
%     set(f,'PaperPositionMode','Auto')
%     maxwindows(f);
%     if savefigures
%         saveas(gcf,['savedfigs/DbFitAbs_' studyID '_' plotID_flow '_plot' num2str(fn) '.fig'],'fig')
%         saveas(gcf,['savedfigs/DbFitAbs_' studyID '_' plotID_flow '_plot' num2str(fn) '.eps'],'epsc2')
%         saveas(gcf,['savedfigs/DbFitAbs_' studyID '_' plotID_flow '_plot' num2str(fn) '.png'],'png')
%     end
% end

%%%%%%%%%%


% tt=['save plotDbFitAbs_' studyID '_' plotID '.mat names avgfit studyIDs analysisIDs sources baselinemarks markstoshow markstolabel excludemarks t0 frameperiod region fMarks BFIleft BFIleft_std BFIright BFIright_std BFIleftparietal BFIleftparietal_std BFIrightparietal BFIrightparietal_std BFIforehead BFIforehead_std BFIparietal BFIparietal_std rBFIforehead rBFIparietal'];
% eval(tt);
% for fn = 1:3
%     subplot(1,3,fn)
%     set(gca,'FontSize',25)
%     set(findobj(gca,'Type','text'),'FontSize',25)
% %     xlabh = get(gca,'XLabel');
% %     set(xlabh,'Position',get(xlabh,'Position') - [0 .2 0]);    
% end
set(f,'PaperPositionMode','Auto')
maxwindows(f);
fn = 1;
if savefigures
    saveas(gcf,['savedfigs/RegAbsProps_' studyID '_' plotID_oxy '_plot' num2str(fn) '.fig'],'fig')
    saveas(gcf,['savedfigs/RegAbsProps_' studyID '_' plotID_oxy '_plot' num2str(fn) '.eps'],'epsc2')
    saveas(gcf,['savedfigs/RegAbsProps_' studyID '_' plotID_oxy '_plot' num2str(fn) '.png'],'png')
end
