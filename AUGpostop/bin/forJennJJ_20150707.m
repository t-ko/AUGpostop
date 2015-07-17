%% Data Comparison with Jenn by Tiffany
% 07-07-2015
% run after PlotDbFitAbs for postop AUG data
close all
plotfigs = 0;
plotID_flow='BFIabs_fitavg';
plotID_oxy='absprops';

clear ALL_mualeft ALL_muaright ALL_muspleft ALL_muspright
% Weight optical prop average by date
for i = 1:length(ALLstudyIDs)
    studyID = ALLstudyIDs{i};        
    tt=[ 'load plotAbsProps_' studyID '_absprops.mat' ];
    eval(tt);
    for RR = 1:length(analysisIDs)
        absID = analysisIDs{RR};
        ext='';
        fname1=[ studyID '_' absID ext '_baselines.mat'];
        load(fname1);
        fdir=['..' filesep studyID filesep absID filesep];
        files=dir([ fdir 'Data_*.txt']);
        fname=files(1).name(6:end-4);
        if length(fname)>13
            DateVector(RR,:) = datevec([fname(5:6) '' filesep fname(7:8) '' filesep fname(1:4) ' ' fname(end-5:end-4) ':' fname(end-3:end-2) ':' fname(end-1:end)]);
        else
            DateVector(RR,:) = datevec([fname(1:2) '' filesep fname(3:4) '' filesep fname(5:6) ' ' fname(end-5:end-4) ':' fname(end-3:end-2) ':' fname(end-1:end)]);
        end
        analysis_date(RR) = datenum([DateVector(RR,1:3) 0 0 0]);
    end
    ALL_mualeftavg{i} = DateWeightedProp(analysis_date, mualeft).*ones(size(mualeft));
    ALL_muarightavg{i} = DateWeightedProp(analysis_date, muaright).*ones(size(muaright));
    ALL_muspleftavg{i} = DateWeightedProp(analysis_date, muspleft).*ones(size(muspleft));
    ALL_musprightavg{i} = DateWeightedProp(analysis_date, muspright).*ones(size(muspright));
    ALL_mualeft{i} = mualeft;
    ALL_muaright{i} = muaright;
    ALL_muspleft{i} = muspleft;
    ALL_muspright{i} = muspright;
    clear mualeft muaright muspleft muspright analysis_date
end


col_names = {'ID','Time','Hb','HbO2','THC','StO2','BFI','musp_forehead','mua_forehead','musp_leftavg','musp_rightavg','mua_leftavg','mua_rightavg','musp_left','musp_right','mua_left','mua_right'};
for i = 1:length(ALLstudyIDs)
    ID = ALLstudyIDs{i};
    ALL_studyID{i} = str2num(ID(4:6)).*ones(1,length(ALL_time_oxy{i}));
    fname2=[ 'plotAbsProps_' ID '_' plotID_oxy '.mat'];
    load([ fname2 ]); 
    t = ALL_time_oxy{i}.*24; %convert days to hours post-op
    if plotfigs
        figure(1); plot(t,ALL_Hbforehead{i},'Color',colors(i,:),'LineWidth',3);
        ylabel(texlabel('[Hb] ({mu}mol/mL)'),'FontSize',25);
        figure(2); plot(t,ALL_HbO2forehead{i},'Color',colors(i,:),'LineWidth',3);
        ylabel(texlabel('[HbO_2] ({mu}mol/mL)'),'FontSize',25);
        figure(3); plot(t,ALL_THCforehead{i},'Color',colors(i,:),'LineWidth',3);
        ylabel(texlabel('Total Hgb Concentration {mu}mol/mL'),'FontSize',25);
        figure(4); plot(t,ALL_StO2forehead{i},'Color',colors(i,:),'LineWidth',3);
        ylabel(texlabel('StO_2 (%)'),'FontSize',25);
        figure(5); plot(t,ALL_BFIforehead{i},'Color',colors(i,:),'LineWidth',3);
        ylabel(texlabel('BFI (cm^2/s)'),'FontSize',25);
        figure(6); plot(t,ALL_muspforehead{i},'Color',colors(i,:),'LineWidth',3);
        ylabel(texlabel('mu_s'' (cm^{-1})'),'FontSize',25);
        figure(7); plot(t,ALL_muaforehead{i},'Color',colors(i,:),'LineWidth',3);
        ylabel(texlabel('mu_a (cm^{-1})'),'FontSize',25);
    end
end

if plotfigs
    for fig = 1:7
        f = figure(fig); 
        xlabel('Post-Operative Day','FontSize',25)
        set(gca,'FontSize',25)
        grid on
        set(f,'PaperPositionMode','Auto')
        maxwindows(f);
    end
end

ID = [ALL_studyID{:}]';
Time = [ALL_time_oxy{:}]';
Hb = [ALL_Hbforehead{:}]';
HbO2 = [ALL_HbO2forehead{:}]';
THC = [ALL_THCforehead{:}]';
StO2 = [ALL_StO2forehead{:}]';
BFI = [ALL_BFIforehead{:}]';
musp = [ALL_muspforehead{:}]';
mua = [ALL_muaforehead{:}]';
muspleftavg = [ALL_muspleftavg{:}]';
musprightavg = [ALL_musprightavg{:}]';
mualeftavg = [ALL_mualeftavg{:}]';
muarightavg = [ALL_muarightavg{:}]';
muspleft = [ALL_muspleft{:}]';
muspright = [ALL_muspright{:}]';
mualeft = [ALL_mualeft{:}]';
muaright = [ALL_muaright{:}]';

% ALLdata_days = table(ID,Time,Hb,HbO2,THC,StO2,BFI,musp,mua,'VariableNames',col_names)
ALLdata_hours = table(ID,Time.*24,Hb,HbO2,THC,StO2,BFI,musp,mua,muspleftavg,musprightavg,mualeftavg,muarightavg,muspleft,muspright,mualeft,muaright,'VariableNames',col_names)
