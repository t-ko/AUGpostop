%Dbfit code template- EB 6/26/09
%see labnote book on why I chose some of these parameters
%Right now only works for 1 s-d separation

clear all
close all

%% FILENAME PARAMETERS
studyID='aug068_parietal';
prefix='BFI_1';
analysisID='BFI1abs_2layer_noCSF_fitavg';
absID='postop_2';
%Extension added to filename, if nothing, just use ext='';
exten={''};

%% DCS PARAMETERS
s=1; %source select, 1 or 2
usedflowdets=[2 3 4 5 6 7 8];%must be consecutive, cannot be ex. 1,3,4 or 2,4.  must change code if this is the case 
fitavg = 1; %average all usedflowdets during correlation fitting?
% regionmarks = {[3 4]}; % leave empty cell {} if no region discrimination
% regionlabels = {'Left Parietal'};
% regionthickness = [0.37];
regionmarks = {[1 2];[3 4];[5 6];[8 9]}; % leave empty cell {} if no region discrimination
regionlabels = {'Left Forehead';'Left Parietal';'Right Forehead';'Right Parietal'};
regionthickness = [0.47 0.37 0.45 0.31];
marksforgamma = [1 2];
collectedo2=1;%Was oxygenation data collected in this study?  If so, will use changes in mua in our Dbfit
syncmark=1;

%% PROBE PARAMETERS
Ns=1;
Nd=8;
r=2.7; % Custom probe 1

%% FIT CONSTANTS
n0=1.4;%index of refraction for tissue
c=2.99792458e10; %speed of light in vacuum,cm/s
vo=c/n0; %speed of light in medium
lambda=7.85e-5;
R=-1.440./n0^2+0.710/n0+0.668+0.0636.*n0;
ze=2./3.*(1+R)./(1-R); %check this number, no need for 1/musp as Cecil takes ca\re of it in the files
muao=0.1;
muspo=7;
Dbr0=1e-8; %background brownian diffusion coeff
D=vo./(3.0.*muspo); %background diffusion coeff
k0=2*pi*n0/lambda; %this is the k0 for flow!

%% FITTING SETTINGS:
fitbeta=1; %=0 do not fit beta,==1 fit beta, DO NOT RECOMMEND FITTING BETA UNLESS INTENSITY > 20kHz
betafitmin=0.4;
betafitmax=0.55;
op=optimset('fminsearch');
options=optimset(op,'MaxIter',300000,'MaxFunEvals',300000,'TolFun',1.000e-16,'TolX',1.000e-16,'Display','Final');
startcorr=6;% if you are fitting beta, want to use a higher number here (e.g. 9), if not fitting beta, use smallest number possible
datalength=90; %chao's way of cutting data for decay.
avgnum=6;%How many points to average in each curve before fitting
cutoff=1.01;%where to cut correlation curve off when fitting
setbeta=0.50; %used when calculating sigma
x0=[1e-8 setbeta];%Initial guess for Db and beta
cutoffmeanerror_perc=10; %mean error suggested cutoff, upper threshold; whole integer percentage removed

%%%%PLOT SETTINGS
intenymax=70;%ylim for intensity plot

%Load flow data
minfiles=0;
measnum=0;
morefiles=1; %1:more files to read, 0:quit
snum=1;

%Step 1, plot amplitudes for each file extension and assess probe off marks
%and baseline marks
measnumtmp=0;
measnumtmp1=0;
for e=1:length(exten)
    fdir1 = ['../' studyID '/'];
    fname1 = [ prefix '_' ];
    if ~isempty(exten{e})
        fname1 = [ fname1 char(exten(e)) '_'];
    end

    %Load flow data
    minfiles=0;
    measnum=0;
    morefiles=1; %1:more files to read, 0:quit
    snum=1;

    %Load all data
    while (morefiles)

        if snum>Ns
            snum=1;
            measnum=measnum+1;
        end

        numfiles=measnum.*Ns+snum+minfiles-1;
        fname=[fdir1 fname1 'flow_' sprintf('%01d',numfiles) '.dat'];

        if exist(fname)==2 %file exists
            %Find time of frame:
            fid = fopen([ fname ], 'r');
            [tmpdata, count]=fscanf(fid,'%c %s',3);
            clear tmpdata
            [tmpdata, count]=fscanf(fid,'%c %s',2);
            time{measnum+1}=tmpdata;
            fclose(fid);
            
            %Load correlator data
            data=load(fname);
            
            %Record marks
            if data(end,1)>0
                Marksflow(data(end,1))=measnum+1;
            end

            if size(data,2)==9
                intensitydata(snum,measnum+1,:)=data(1,2:9);
                corrs(snum,measnum+1,:,:)=data(3:end,2:9);
            else
                intensitydata(snum,measnum+1,:)=data(1,2:5);
                corrs(snum,measnum+1,:,:)=data(3:end,2:5);
            end
            %Find frames with light leakage
            for d=1:size(corrs,4)
                if mean(corrs(snum,measnum+1,datalength:datalength+20,d))>1.01
                    corrs(snum,measnum+1,:,d)=NaN;
                end
            end
            %Find Beta
            Betasave(measnum+1,:)=nanmean(squeeze(corrs(snum,measnum+1,usedflowdets)),2)-1;
            
            taus=data(3:end,1);
            snum=snum+1;
        elseif exist(fname)==0 %file does not exist
            morefiles=0;
        end
    end
    marksperext(e)=length(Marksflow);
    filelength(e)=size(intensitydata,2)+measnumtmp;
    %Make running amplitude and Marks variables
    intensitydatatmp(:,measnumtmp+1:filelength(e),:)=intensitydata;
    corrstmp(:,measnumtmp+1:filelength(e),:,:)=corrs;
    timetmp(:,measnumtmp+1:filelength(e))=time;
    Marksflowtmp(measnumtmp1+1:length(Marksflow)+measnumtmp1)=Marksflow+measnumtmp+1;
    measnumtmp=filelength(e);
    measnumtmp1=length(Marksflow)+measnumtmp1;
    %Get integration time (sec)
    t=data(1,1)/1000;

    clear data Marksflow corrs intensitydata time
end
corrstmp0=corrstmp; %placeholder storage variable
Marksflow=Marksflowtmp;
if Marksflow(end) >= size(corrstmp0,2)
    Marksflow(end) = size(corrstmp0,2)-1;
end

%Calculate time points corresponding to each flow curve
timeaxis_flow=datenum(timetmp,'HH:MM:SS')-floor(datenum(timetmp,'HH:MM:SS'));%In arbitrary units--1a.u.=24hrs, counting from 1/1/2000.

%Data in DCS files is not in military time, so add 12 hours to all data if
%need be
if timeaxis_flow(1)<0.5
    timeaxis_flow=timeaxis_flow+0.5;
end

%Unwrap time vector
for i=1:length(timeaxis_flow)-1
    if abs(timeaxis_flow(i+1)-timeaxis_flow(i))>0.4
        timeaxis_flow(i+1)=timeaxis_flow(i+1)+0.5;
    end
end

%Zero time axis so that t=0 is first flow curve
timeaxis_flow=(timeaxis_flow-timeaxis_flow(1))*1.44e3;%Convert time to minutes

%Correct for time since off bypass
timeaxis_flow=timeaxis_flow-timeaxis_flow(Marksflow(syncmark));
figure,plot(timeaxis_flow)

%Determine bin width for each tau
T=zeros(size(taus));
for indt=1:length(T)-1
    T(indt)=taus(indt+1)-taus(indt);
end


% if length(usedflowdets)==1
%     fitavg=0;
% else
%     %Ask for user input as to whether to fit average or individual curves
%     fitavg=input('Fit average of corr. curves? 0=no, 1=yes  ');
% end

    %% Import each detector data separately
    intensitydata=squeeze(intensitydatatmp(s,:,:));
    corrs=squeeze(corrstmp0(s,:,:,:));

    figure,plot(intensitydata(:,usedflowdets),'.-','MarkerSize',20,'LineWidth',3)
    %end
    xlabel('Frame','FontSize',25)
    ylabel('Intensity (kHz)','FontSize',25)
    axis tight
    ylim([0 intenymax])
    grid on
    set(gca,'FontSize',20)
    set(gcf,'PaperPositionMode','Auto')
    saveas(gcf,['../' studyID '/savedfigs/Intensity_' studyID '_S' num2str(s) '_' analysisID '.fig'],'fig')
    saveas(gcf,['../' studyID '/savedfigs/Intensity_' studyID '_S' num2str(s) '_' analysisID '.jpg'],'jpg')

    %% Import absolute mua and musp, if measured (TK 2014-8-2)
    mua=muao+zeros(1,size(intensitydata,1));
    musp=muspo+zeros(1,size(intensitydata,1));;%select scattering at wavelength closest to 785 (TK 2014-08)
    if exist([ studyID '_' absID '_baselines.mat'],'file')
        tt=['load ' studyID '_' absID '_baselines.mat'];
        eval(tt);
        for rn = 1:length(regionmarks)
            rmarks = Marksflow(regionmarks{rn});
            if ~isempty(findstr(regionlabels{rn},'eft')) && ~isempty(findstr(regionlabels{rn},'orehead'))
                mua(rmarks(1):rmarks(2)) = mean(leftmua(:,3));
                musp(rmarks(1):rmarks(2)) = mean(leftmusp(:,3));
            elseif ~isempty(findstr(regionlabels{rn},'eft')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
                mua(rmarks(1):rmarks(2)) = mean(leftparietalmua(:,3));
                musp(rmarks(1):rmarks(2)) = mean(leftparietalmusp(:,3));
            elseif ~isempty(findstr(regionlabels{rn},'ight')) && ~isempty(findstr(regionlabels{rn},'orehead'))
                mua(rmarks(1):rmarks(2)) = mean(rightmua(:,3));
                musp(rmarks(1):rmarks(2)) = mean(rightmusp(:,3));
            elseif ~isempty(findstr(regionlabels{rn},'ight')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
                mua(rmarks(1):rmarks(2)) = mean(rightparietalmua(:,3));
                musp(rmarks(1):rmarks(2)) = mean(rightparietalmusp(:,3));
            end
        end
        
    end
    
    fitmultframes=input('Fit average of multiple frames?  0=no, 1=yes (this option is for very noisy data)  ');
    if fitmultframes==1
        numframestoavg=input('How many frames would you like to average?   ');
    else
        numframestoavg=1;
    end

    %Define a cutoff intensity to get rid of bad frames
    cutoffintensity=input('Cutoff intensity value (kHz)? (All frames with intensity above this value will be discarded)  ');

    %Plot all correlation curves to get a sense of what beta should be
    if length(usedflowdets)==1
        d=usedflowdets;
        %Find frames where I is bigger than cutoff
        int=find(intensitydata(:,d)>cutoffintensity);
        %Set frames with I above cutoff = NaN
        corrs(int,:,d)=NaN;

        figure,semilogx(taus,corrs(1,:,d),'r-','LineWidth',1)
        for m=2:size(corrs,1)
            hold on,semilogx(taus,corrs(m,:,d),'r-','LineWidth',1)
        end
        hold on,semilogx(taus,nanmean(corrs(:,:,d),1),'k-','LineWidth',5)
        xlabel('\tau','FontSize',25)
        ylabel('g2','FontSize',25)
        title(['Detector ' num2str(d) ],'FontSize',25)
        ylim([0.9 1.6])
        xlim([0 1e-2])
        set(gca,'FontSize',20)
        h=line([taus(datalength) taus(datalength)],[0.9 1.6]);
        set(h,'Color',[0 0 0]);
        set(gcf,'PaperPositionMode','Auto')
        saveas(gcf,['../' studyID '/savedfigs/corrs_det' num2str(d) '_' studyID '_S' num2str(s) '_' analysisID '.fig'],'fig')
        saveas(gcf,['../' studyID '/savedfigs/corrs_det' num2str(d) '_' studyID '_S' num2str(s) '_' analysisID '.jpg'],'jpg')
        clear int
    else
        for d=usedflowdets
            %Find frames where I is bigger than cutoff
            int=find(intensitydata(:,d)>cutoffintensity);
            %Set frames with I above cutoff = NaN
            corrs(int,:,d)=NaN;

            figure,semilogx(taus,corrs(1,:,d),'r-','LineWidth',1)
            for m=2:size(corrs,1)
                hold on,semilogx(taus,corrs(m,:,d),'r-','LineWidth',1)
            end
            hold on,semilogx(taus,nanmean(corrs(:,:,d),1),'k-','LineWidth',5)
            xlabel('\tau','FontSize',25)
            ylabel('g2','FontSize',25)
            title(['Detector ' num2str(d) ],'FontSize',25)
            ylim([0.9 1.6])
            xlim([0 1e-2])
            set(gca,'FontSize',20)
            h=line([taus(datalength) taus(datalength)],[0.9 1.6]);
            set(h,'Color',[0 0 0]);
            set(gcf,'PaperPositionMode','Auto')
            saveas(gcf,['../' studyID '/savedfigs/corrs_det' num2str(d) '_' studyID '_S' num2str(s) '_' analysisID '.fig'],'fig')
            saveas(gcf,['../' studyID '/savedfigs/corrs_det' num2str(d) '_' studyID '_S' num2str(s) '_' analysisID '.jpg'],'jpg')
            clear int
        end
    end

    %Calculate gamma for sigma calculation
    %Will use the mean of all correlation curves taken over the range set above.  In
    %this way, the curve will be smooth, and it will provide us with
    %approximately the right order of magnitude for gamma
    if length(usedflowdets==1)
        corrsavgtmp=nanmean(corrs(:,:,usedflowdets),3);%First take avg over all dets
    else
        corrsavgtmp=nanmean(corrs(:,:,usedflowdets),3);%First take avg over all dets
    end

    diff=abs(nanmean(corrsavgtmp(Marksflow(marksforgamma(1)):Marksflow(marksforgamma(2)),:),1)-(1+setbeta*1/exp(1)));
    ind=find(diff==min(diff));
    gamma=1/taus(min(ind));%use min(ind) in case ind is not a 1x1 vector.
