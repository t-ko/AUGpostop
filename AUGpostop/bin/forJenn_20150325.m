%% Data Comparison with Jenn by Tiffany
% 03-25-2015
% run after PlotDbFitAbs for postop AUG data
close all

col_names = {'ID','Time','Hb','HbO2','THC','StO2','musp','BFI (1e-9)','mua'};
for i = 1:length(IDs)
    ALL_studyID{i} = IDs(i).*ones(1,length(ALL_time_oxy{i}));
    t = ALL_time_oxy{i}.*24; %convert days to hours post-op
    figure(1); plot(t,ALL_Hbforehead{i},'Color',colors(i,:),'LineWidth',3);
    ylabel(texlabel('[Hb] ({mu}mol/mL)'),'FontSize',25);
    figure(2); plot(t,ALL_HbO2forehead{i},'Color',colors(i,:),'LineWidth',3);
    ylabel(texlabel('[HbO_2] ({mu}mol/mL)'),'FontSize',25);
    figure(3); plot(t,ALL_THCforehead{i},'Color',colors(i,:),'LineWidth',3);
    ylabel(texlabel('Total Hgb Concentration {mu}mol/mL'),'FontSize',25);
    figure(4); plot(t,ALL_StO2forehead{i},'Color',colors(i,:),'LineWidth',3);
    ylabel(texlabel('StO_2 (%)'),'FontSize',25);
    figure(5); plot(t,ALL_muspforehead{i},'Color',colors(i,:),'LineWidth',3);
    ylabel(texlabel('mu_s'' (cm^{-1})'),'FontSize',25);
    figure(6); plot(t,ALL_BFIforehead{i},'Color',colors(i,:),'LineWidth',3);
    ylabel(texlabel('BFI (cm^2/s)'),'FontSize',25);
    figure(7); plot(t,ALL_muaforehead{i},'Color',colors(i,:),'LineWidth',3);
    ylabel(texlabel('mu_a (cm^{-1})'),'FontSize',25);
end
for fig = 1:7
    f = figure(fig); 
    xlabel('Post-Operative Day','FontSize',25)
    set(gca,'FontSize',25)
    grid on
    set(f,'PaperPositionMode','Auto')
    maxwindows(f);
end

ID = [ALL_studyID{:}]';
Time = [ALL_time_oxy{:}]';
Hb = [ALL_Hbforehead{:}]';
HbO2 = [ALL_HbO2forehead{:}]';
THC = [ALL_THCforehead{:}]';
StO2 = [ALL_StO2forehead{:}]';
musp = [ALL_muspforehead{:}]';
BFI = [ALL_BFIforehead{:}]';
mua = [ALL_muaforehead{:}]';

ALLdata = table([ID,Time,Hb,HbO2,THC,StO2,musp,BFI,mua],'VariableNames',col_names)
