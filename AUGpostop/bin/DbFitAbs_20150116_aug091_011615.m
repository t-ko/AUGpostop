%Dbfit code template- EB 6/26/09
%see labnote book on why I chose some of these parameters
%Right now only works for 1 s-d separation

clear all
close all

%% FILENAME PARAMETERS
studyID='AUG091_011615';
prefix='BFI_7';
analysisID='BFI_011615abs_fitavg';
absID='absolutes_011615';
%Extension added to filename, if nothing, just use ext='';
exten={''};

%% DCS PARAMETERS
s=1; %source select, 1 or 2
usedflowdets=[6 7 8];%must be consecutive, cannot be ex. 1,3,4 or 2,4.  must change code if this is the case 
fitavg = 1; %average all usedflowdets during correlation fitting?
regionmarks = {[1 2],[3 4]}; % leave empty cell {} if no region discrimination
regionlabels = {'Left Forehead','Right Forehead'};
marksforgamma = [1 2];
collectedo2=1;%Was oxygenation data collected in this study?  If so, will use changes in mua in our Dbfit
syncmark=1;

%% PROBE PARAMETERS
% Custom Probe 1
Ns=1;
Nd=8;
SDS=[2.75 2.75 2.75 2.75 2.75 2.75 2.75 2.75];

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
betafitmin=0.40;
betafitmax=0.55;
op=optimset('fminsearch');
options=optimset(op,'MaxIter',300000,'MaxFunEvals',300000,'TolFun',1.000e-16,'TolX',1.000e-16,'Display','Final');
startcorr=6;% if you are fitting beta, want to use a higher number here (e.g. 9), if not fitting beta, use smallest number possible
datalength=90; %chao's way of cutting data for decay.
avgnum=6;%How many points to average in each curve before fitting
cutoff=1.2;%where to cut correlation curve off when fitting
setbeta=0.50; %used when calculating sigma
x0=[1e-8 setbeta];%Initial guess for Db and beta
cutoffmeanerror_perc=10; %mean error suggested cutoff, upper threshold; whole integer percentage removed

%%%%PLOT SETTINGS
intenymax=100;%ylim for intensity plot

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
    fdir1 = ['..' filesep studyID filesep];
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
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'Intensity_' studyID '_S' num2str(s) '_' analysisID '.fig'],'fig')
    saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'Intensity_' studyID '_S' num2str(s) '_' analysisID '.jpg'],'jpg')

    mua=muao+zeros(1,size(intensitydata,1));
    musp=muspo+zeros(1,size(intensitydata,1));;%select scattering at wavelength closest to 785 (TK 2014-08)
    %% Import patient average mua/musp (TK 2015-05-07)
    if exist([ 'plotAbsProps_' studyID '_absprops.mat' ],'file')
        tt=[ 'load plotAbsProps_' studyID '_absprops.mat' ];
        eval(tt);

        % Weight optical prop average by date
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
        % Search regionlabels for region information, forehead/parietal
        for rn = 1:length(regionmarks)
            rmarks = Marksflow(regionmarks{rn});
            if ~isempty(findstr(regionlabels{rn},'eft')) && ~isempty(findstr(regionlabels{rn},'orehead'))
                mua(rmarks(1):rmarks(2)) = DateWeightedProp(analysis_date, mualeft);
                musp(rmarks(1):rmarks(2)) = DateWeightedProp(analysis_date, muspleft);
                sprintf('Left Forehead - mua: %s, musp: %s', num2str(mua(rmarks(1))), num2str(musp(rmarks(1))))
            elseif ~isempty(findstr(regionlabels{rn},'eft')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
                mua(rmarks(1):rmarks(2)) = DateWeightedProp(analysis_date, mualeftparietal);
                musp(rmarks(1):rmarks(2)) = DateWeightedProp(analysis_date, muspleftparietal);
                sprintf('Left Parietal - mua: %s, musp: %s', num2str(mua(rmarks(1))), num2str(musp(rmarks(1))))
            elseif ~isempty(findstr(regionlabels{rn},'ight')) && ~isempty(findstr(regionlabels{rn},'orehead'))
                mua(rmarks(1):rmarks(2)) = DateWeightedProp(analysis_date, muaright);
                musp(rmarks(1):rmarks(2)) = DateWeightedProp(analysis_date, muspright);
                sprintf('Right Forehead - mua: %s, musp: %s', num2str(mua(rmarks(1))), num2str(musp(rmarks(1))))
            elseif ~isempty(findstr(regionlabels{rn},'ight')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
                mua(rmarks(1):rmarks(2)) = DateWeightedProp(analysis_date, muarightparietal);
                musp(rmarks(1):rmarks(2)) = DateWeightedProp(analysis_date, musprightparietal);
                sprintf('Right Parietal - mua: %s, musp: %s', num2str(mua(rmarks(1))), num2str(musp(rmarks(1))))
            end
        end 
    elseif exist([ studyID '_' absID '_baselines.mat'],'file')
    %% Import absolute mua and musp, if measured (TK 2014-8-2)
        tt=['load ' studyID '_' absID '_baselines.mat'];
        eval(tt);
        for rn = 1:length(regionmarks)
            rmarks = Marksflow(regionmarks{rn});
            if ~isempty(findstr(regionlabels{rn},'eft')) && ~isempty(findstr(regionlabels{rn},'orehead'))
                mua(rmarks(1):rmarks(2)) = nanmean(leftmua(:,3));
                musp(rmarks(1):rmarks(2)) = nanmean(leftmusp(:,3));
            elseif ~isempty(findstr(regionlabels{rn},'eft')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
                mua(rmarks(1):rmarks(2)) = nanmean(leftparietalmua(:,3));
                musp(rmarks(1):rmarks(2)) = nanmean(leftparietalmusp(:,3));
            elseif ~isempty(findstr(regionlabels{rn},'ight')) && ~isempty(findstr(regionlabels{rn},'orehead'))
                mua(rmarks(1):rmarks(2)) = nanmean(rightmua(:,3));
                musp(rmarks(1):rmarks(2)) = nanmean(rightmusp(:,3));
            elseif ~isempty(findstr(regionlabels{rn},'ight')) && (~isempty(findstr(regionlabels{rn},'ietal')) || ~isempty(findstr(regionlabels{rn},'ide')))
                mua(rmarks(1):rmarks(2)) = nanmean(rightparietalmua(:,3));
                musp(rmarks(1):rmarks(2)) = nanmean(rightparietalmusp(:,3));
            end
        end        
    end
    
    % Define frames with good data to calculate intensity cutoff
    ONframes = [];
    for ii = 1:length(regionmarks)
        try
            ONframes = [ONframes Marksflow(regionmarks{ii}(1)):(Marksflow(regionmarks{ii}(2)))];
        catch
        end
    end
    %Define a cutoff intensity to get rid of bad frames
%     cutoffintensity=input('Cutoff intensity value (kHz)? (All frames with intensity above this value will be discarded)  ');
    IQR = quantile(intensitydata(ONframes,usedflowdets),3);    
    try 
        p75_intensity = max(IQR(3,:))
        cutoffintensity = max((IQR(3,:)-IQR(1,:))+IQR(3,:)) %throw away outlier intensity
    catch
        p75_intensity = IQR(3)
        cutoffintensity = (IQR(3)-IQR(1))+IQR(3) %throw away outlier intensity
    end
    %     fitmultframes=input('Fit average of multiple frames?  0=no, 1=yes (this option is for very noisy data)  ');
    numframestoavg = max([min([round(100/p75_intensity) 3]) 1])
    if numframestoavg > 1
        fitmultframes=1;
    else
        fitmultframes=0;
    end
%     if fitmultframes==1
%         numframestoavg=input('How many frames would you like to average?   ');
%     else
%         numframestoavg=1;
%     end

    % Redefine ONframes to account for frame averaging
    ONframes = [];
    for ii = 1:length(regionmarks)
        try
            ONframes = [ONframes Marksflow(regionmarks{ii}(1)):(Marksflow(regionmarks{ii}(2))-numframestoavg)];
        catch
        end
    end
    % Remove data above intensity cutoff
    d=usedflowdets;
    %Find frames where I is bigger than cutoff
    int=find(intensitydata(:,d)>cutoffintensity);
    %Set frames with I above cutoff = NaN
    corrs(int,:,d)=NaN;
    clear int
    
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

    %Now fit each frame of data
    figure()

    for RR = 1:length(regionlabels)
        numframestoavg = Marksflow(regionmarks{RR}(2))-Marksflow(regionmarks{RR}(1));
        if fitavg==1
            r = unique(SDS(s,usedflowdets));
            
            %% OPTIMIZE BETA
            measnum = Marksflow(regionmarks{RR}(1))-1;            
            %First take avg over all dets
            corrsavgtmp=nanmean(corrs(:,:,usedflowdets),3);%First take avg over all dets
            
            % PLOT G2
            figure,semilogx(taus,nanmean(corrsavgtmp(measnum+1:measnum+numframestoavg,:),1),'k-','LineWidth',5)
            for m=measnum+1:measnum+numframestoavg
                hold on,semilogx(taus,corrsavgtmp(m,:),'r-','LineWidth',1)
            end
            xlabel('\tau','FontSize',25)
            ylabel('g2','FontSize',25)
            title(['Average g2: ' regionlabels{RR} ],'FontSize',25)
            ylim([0.9 1.6])
            xlim([0 1e-2])
            set(gca,'FontSize',20)
            h=line([taus(datalength) taus(datalength)],[0.9 1.6]);
            set(h,'Color',[0 0 0]);
            set(gcf,'PaperPositionMode','Auto')
            
            %% Calculate gamma for sigma calculation
            %Will use the mean of all correlation curves taken over the range set above.  In
            %this way, the curve will be smooth, and it will provide us with
            %approximately the right order of magnitude for gamma
            if length(usedflowdets==1)
                corrsavgtmp=nanmean(corrs(:,:,usedflowdets),3);%First take avg over all dets
            else
                corrsavgtmp=nanmean(corrs(:,:,usedflowdets),3);%First take avg over all dets
            end
            corrsdiff=abs(nanmean(corrsavgtmp(Marksflow(marksforgamma(1)):Marksflow(marksforgamma(2)),:),1)-(1+setbeta*1/exp(1)));
            ind=find(corrsdiff==min(corrsdiff));
            gamma=1/taus(min(ind));%use min(ind) in case ind is not a 1x1 vector.

            intensityavg(measnum+1)=nanmean(nanmean(intensitydata(measnum+1:measnum+numframestoavg,usedflowdets),2),1);
            corrsavg(measnum+1,:)=nanmean(nanmean(corrs(measnum+1:measnum+numframestoavg,:,usedflowdets),3),1);
            
            foo=min(find(slidingavg(corrsavg(measnum+1,:),avgnum)<=cutoff));
            if isempty(foo) || foo<startcorr
                tmpf=datalength;
            else
                tmpf=min(find(slidingavg(corrsavg(measnum+1,:),avgnum)<=cutoff));
            end
            corrstmpavg=squeeze(corrsavg(measnum+1,startcorr:tmpf));
            corrstmpavg=slidingavg(corrstmpavg,avgnum);
            taustmp=taus(startcorr:tmpf);
            %Calculate noise from Chao's noise model
            sigma=1./intensityavg(measnum+1).*sqrt(1./T./t).*sqrt(1+mean(setbeta).*exp(-gamma*taus));
            %Fit both Beta and BFI at the same time, bound Beta fit to
            %reasonable values
            [betaDbfitavg(:,measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_betaandDB_new1_withsigma,x0,...
            [0 betafitmin],[1e-1 betafitmax],options,r,taustmp,musp(measnum+1),mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);
            regionbeta(RR)=squeeze(betaDbfitavg(2,measnum+1));
            setbeta = regionbeta(RR)
            text(5e-7,1,['I=' num2str(intensityavg(measnum+1),'%6.2f') ],'FontSize',16 )
            text(5e-7,1.1,['\beta=' num2str(setbeta,'%6.2f') ],'FontSize',16 )
            text(5e-7,1.2,['Db=' num2str(squeeze(betaDbfitavg(1,measnum+1)),'%6.2e') ],'FontSize',16 )             
            
            saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'corrsavg_' studyID '_S' num2str(s) '_' analysisID '_' regionlabels{RR} '.fig'],'fig')
            saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'corrsavg_' studyID '_S' num2str(s) '_' analysisID '_' regionlabels{RR} '.jpg'],'jpg')

            
            numframestoavg=1;
            
            %% FIX BETA FITTING        
            for measnum=0:size(intensitydata,1)-numframestoavg %UPDATED (TK) 02-22-15 to only include marked frames
                if (measnum+numframestoavg >= Marksflow(regionmarks{RR}(1))) && (measnum+numframestoavg <= Marksflow(regionmarks{RR}(2)))
                    intensityavg(measnum+1)=nanmean(nanmean(intensitydata(measnum+1:measnum+numframestoavg,usedflowdets),2),1);
                    corrsavg(measnum+1,:)=nanmean(nanmean(corrs(measnum+1:measnum+numframestoavg,:,usedflowdets),3),1);
                    if isempty(find(isnan(corrsavg(measnum+1,:)))) & ~isnan(mua(measnum+1))
                        foo=min(find(slidingavg(corrsavg(measnum+1,:),avgnum)<=cutoff));
                        if isempty(foo) || foo<startcorr
                            tmpf=datalength;
                        else
                            tmpf=min(find(slidingavg(corrsavg(measnum+1,:),avgnum)<=cutoff));
                        end
                        corrstmpavg=squeeze(corrsavg(measnum+1,startcorr:tmpf));
                        corrstmpavg=slidingavg(corrstmpavg,avgnum);
                        taustmp=taus(startcorr:tmpf);
                        %Calculate noise from Chao's noise model
                        sigma=1./intensityavg(measnum+1).*sqrt(1./T./t).*sqrt(1+mean(setbeta).*exp(-gamma*taus));

                        %Dont fit beta, bound Db fit to reasonable values
                        Betasaveavg(measnum+1)=mean(setbeta);
                        [Dbfitavg(measnum+1),fvalavg(measnum+1),exitflagavg(measnum+1)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
                        Betasaveavg(measnum+1),r,taustmp,musp(measnum+1),mua(measnum+1),k0,ze,corrstmpavg,length(taustmp),sigma);

                        g1avg(:,measnum+1)=sqrt(abs((corrsavg(measnum+1,:)-1)./Betasaveavg(measnum+1)));
                        Curvefitavg(:,measnum+1)=g1fitx(Dbfitavg(measnum+1),r,taus,musp(measnum+1),mua(measnum+1),k0,ze);
                        Curvefitg2avg(:,measnum+1)=g2fitx([Dbfitavg(measnum+1) Betasaveavg(measnum+1)],r,taus,musp(measnum+1),mua(measnum+1),k0,ze);
                        
                        %Calculate error in fit
                        indtmp(measnum+1)=min(find(abs(squeeze(Curvefitavg(:,measnum+1))-0.3)==min(abs(squeeze(Curvefitavg(:,measnum+1))-0.3))));%Use min in case size(ind)>1
                        errorfitavg(:,measnum+1)=(g1avg(:,measnum+1)-squeeze(Curvefitavg(:,measnum+1)))./squeeze(Curvefitavg(:,measnum+1))*100;
                        meanerroravg(measnum+1)=mean(errorfitavg(1:indtmp(measnum+1),measnum+1));
                        stderroravg(measnum+1)=std(errorfitavg(1:indtmp(measnum+1),measnum+1));
                        
                        %% Plot g1 fit
                        figure(gcf);
                        hold off
                        semilogx(taus(1:size(g1avg,1)),g1avg(:,measnum+1))
                        hold on; axis tight;
                        semilogx(taus,Curvefitavg(:,measnum+1),'black')
                        title({['g1 Semi-Infinite AvgFit: Frame ' num2str(measnum+1)],['mua: ' num2str(mua(measnum+1)) ', musp: ' num2str(musp(measnum+1))]},'FontSize', 16);
                        text(5e-7,0.1,['I=' num2str(intensityavg(measnum+1),'%6.2f') ],'FontSize',16 )
                        text(5e-7,0.2,['\beta=' num2str(Betasaveavg(measnum+1),'%6.2f') ],'FontSize',16 )
                        xlim([0 10^-3])
                        ylim([0 1.1])
                        set(gca, 'FontSize', 16)
                        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'fit_' studyID '_S' num2str(s) '_' analysisID '_f' num2str(measnum+1) '.fig'],'fig')
                        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'fit_' studyID '_S' num2str(s) '_' analysisID '_f' num2str(measnum+1) '.jpg'],'jpg')

                    else
                        indtmp(measnum+1)=NaN;
                        errorfitavg(:,measnum+1)=ones(1,length(taus)).*NaN;
                        meanerroravg(measnum+1)=NaN;
                        stderroravg(measnum+1)=NaN;
                        Betasaveavg(measnum+1)=NaN;
                        Dbfitavg(measnum+1)=NaN;
                        fvalavg(measnum+1)=NaN;
                        exitflagavg(measnum+1)=NaN;
                        g1avg(:,measnum+1)=ones(1,length(taus)).*NaN;
                        Curvefitavg(:,measnum+1)=ones(1,length(taus)).*NaN;
                        Curvefitg2avg(:,measnum+1)=ones(1,length(taus)).*NaN;
                    end
                    clear corrstmpavg taustmp
                end
            end
        else
            %Fit data from each detector
            for d=1:length(usedflowdets)
                i=usedflowdets(d);
                r=SDS(s,i);
                
                %% Calculate gamma for sigma calculation
                corrsdiff=abs(nanmean(corrs(Marksflow(marksforgamma(1)):Marksflow(marksforgamma(2)),:,i),1)-(1+setbeta*1/exp(1)));
                ind=find(corrsdiff==min(corrsdiff));
                gamma=1/taus(min(ind));%use min(ind) in case ind is not a 1x1 vector.
                
                % OPTIMIZE BETA
                measnum = Marksflow(regionmarks{RR}(1));   
                corrsmean=nanmean(corrs(measnum+1:measnum+numframestoavg,:,i),1);
                                
                % PLOT G2
                figure,semilogx(taus,corrsmean,'k-','LineWidth',5)
                for m=measnum+1:measnum+numframestoavg
                    hold on,semilogx(taus,corrs(m,:,i),'r-','LineWidth',1)
                end
                xlabel('\tau','FontSize',25)
                ylabel('g2','FontSize',25)
                title(['Detector ' num2str(i) ': ' regionlabels{RR} ],'FontSize',25)
                ylim([0.9 1.6])
                xlim([0 1e-2])
                set(gca,'FontSize',20)
                h=line([taus(datalength) taus(datalength)],[0.9 1.6]);
                set(h,'Color',[0 0 0]);
                set(gcf,'PaperPositionMode','Auto')
                
                
                foo=min(find(slidingavg(corrsmean,avgnum)<=cutoff));
                if isempty(foo) || foo<startcorr
                    tmpf(i)=datalength;
                else
                    tmpf(i)=min(find(slidingavg(corrsmean,avgnum)<=cutoff));
                end
                corrstmp=squeeze(corrsmean(startcorr:tmpf(i)));
                corrstmp=slidingavg(corrstmp,avgnum);
                taustmp=taus(startcorr:tmpf(i));
                %Calculate noise from Chao's noise model
                sigma(:,i)=1./intensitydata(measnum+1,i).*sqrt(1./T./t).*sqrt(1+mean(setbeta).*exp(-gamma*taus));
                %Fit both Beta and BFI at the same time, bound Beta fit to
                %reasonable values
                [betaDbfit(:,measnum+1,i),fval(measnum+1,i),exitflag(measnum+1,i)]=fminsearchbnd(@xg2fitx_betaandDB_new1_withsigma,x0,...
                [0 betafitmin],[1e-1 betafitmax],options,r,taustmp,musp(measnum+1),mua(measnum+1),k0,ze,corrstmp,length(corrstmp),sigma(:,i));
                regionbeta(RR,i)=squeeze(betaDbfit(2,measnum+1,i));
                setbeta = regionbeta(RR,i)
                
                text(5e-7,1,['I=' num2str(intensitydata(measnum+1,i),'%6.2f') ],'FontSize',16 )
                text(5e-7,1.1,['\beta=' num2str(setbeta,'%6.2f') ],'FontSize',16 )
                text(5e-7,1.2,['Db=' num2str(squeeze(betaDbfit(1,measnum+1,i)),'%6.2e') ],'FontSize',16 )             
            
                saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'corrs_' studyID '_S' num2str(s) '_' analysisID '_' regionlabels{RR} '_det' num2str(i) '.fig'],'fig')
                saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'corrs_' studyID '_S' num2str(s) '_' analysisID '_' regionlabels{RR} '_det' num2str(i) '.jpg'],'jpg')            
                
                numframestoavg = 1;
                  
                %% FIX BETA FITTING        
                for measnum=0:size(intensitydata,1)-numframestoavg %UPDATED (TK) 02-22-15 to only include marked frames
                    if measnum+numframestoavg >= Marksflow(regionmarks{RR}(1)) && measnum+numframestoavg <= Marksflow(regionmarks{RR}(2)) && isempty(find(isnan(corrs(measnum+1,:,i)))) && ~isnan(mua(measnum+1))
                        corrsmean=nanmean(corrs(measnum+1:measnum+numframestoavg,:,i),1);
                        foo=min(find(slidingavg(corrsmean,avgnum)<=cutoff));
                        if isempty(foo) || foo<startcorr
                            tmpf(i)=datalength;
                        else
                            tmpf(i)=min(find(slidingavg(corrsmean,avgnum)<=cutoff));
                        end
                        corrstmp=squeeze(corrsmean(startcorr:tmpf(i)));
                        corrstmp=slidingavg(corrstmp,avgnum);
                        taustmp=taus(startcorr:tmpf(i));
                        %Calculate noise from Chao's noise model
                        sigma(:,i)=1./intensitydata(measnum+1,i).*sqrt(1./T./t).*sqrt(1+mean(setbeta).*exp(-gamma*taus));
                        %Dont fit beta, bound Db fit to reasonable values
                        Betasave(measnum+1,i)=mean(setbeta);
                        [Dbfit(measnum+1,i),fval(measnum+1,i),exitflag(measnum+1,i)]=fminsearchbnd(@xg2fitx_new1_withsigma,x0(1),[0],[1e-2],options,...
                        Betasave(measnum+1,i),r,taustmp,musp(measnum+1),mua(measnum+1),k0,ze,corrstmp,length(taustmp),sigma(:,i));

                        g1(:,i,measnum+1)=sqrt(abs((corrs(measnum+1:measnum+numframestoavg,:,i)-1)./Betasave(measnum+1,i)));
                        Curvefit(:,measnum+1,i)=g1fitx(Dbfit(measnum+1,i),r,taus,musp(measnum+1),mua(measnum+1),k0,ze);
                        Curvefitg2(:,measnum+1,i)=g2fitx([Dbfit(measnum+1,i) Betasave(measnum+1,i)],r,taus,musp(measnum+1),mua(measnum+1),k0,ze);
                        %Calculate error in fit
                        indtmp(measnum+1,i)=min(find(abs(squeeze(Curvefit(:,measnum+1,i))-0.3)==min(abs(squeeze(Curvefit(:,measnum+1,i))-0.3))));%Use min in case size(ind)>1
                        errorfit(:,measnum+1,i)=(g1(:,i,measnum+1)-squeeze(Curvefit(:,measnum+1,i)))./squeeze(Curvefit(:,measnum+1,i))*100;
                        meanerror(measnum+1,i)=mean(errorfit(1:indtmp(measnum+1,i),measnum+1,i));
                        stderror(measnum+1,i)=std(errorfit(1:indtmp(measnum+1,i),measnum+1,i));
                        %% Plot g1 fit
                        figure(gcf);
                        hold off
                        semilogx(taus(1:size(g1,1)),g1(:,i,measnum+1))
                        hold on; axis tight;
                        semilogx(taus,Curvefit(:,measnum+1,i),'black')
                        title({['g1 Semi-Infinite: Frame ' num2str(measnum+1)],['mua: ' num2str(mua(measnum+1)) ', musp: ' num2str(musp(measnum+1))]},'FontSize', 16);
                        text(5e-7,0.1,['I=' num2str(intensitydata(measnum+1,i),'%6.2f') ],'FontSize',16 )
                        text(5e-7,0.2,['\beta=' num2str(Betasave(measnum+1,i),'%6.2f') ],'FontSize',16 )
                        xlim([0 10^-3])
                        ylim([0 1.1])
                        set(gca, 'FontSize', 16)
                        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'fit_' studyID '_S' num2str(s) '_' analysisID '_det' num2str(i) '_f' num2str(measnum+1) '.fig'],'fig')
                        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'fit_' studyID '_S' num2str(s) '_' analysisID '_det' num2str(i) '_f' num2str(measnum+1) '.jpg'],'jpg')


                    else
                        indtmp(measnum+1,i)=NaN;
                        errorfit(:,measnum+1,i)=ones(1,length(taus)).*NaN;
                        meanerror(measnum+1,i)=NaN;
                        stderror(measnum+1,i)=NaN;
                        g1(:,i,measnum+1)=ones(1,length(taus)).*NaN;
                        Curvefit(:,measnum+1,i)=ones(1,length(taus)).*NaN;
                        Curvefitg2(:,measnum+1,i)=ones(1,length(taus)).*NaN;
                        Betasave(measnum+1,i)=NaN;
                        Dbfit(measnum+1,i)=NaN;
                    end
                    clear corrstmp taustmp
                end
            end
        end
    end
        
                
if ~fitavg
    for d=1:length(usedflowdets)
        i=usedflowdets(d);
        r=SDS(s,i);
        %Make sure fit converged
        ind=find(exitflag(:,i)==0);
        Dbfit(ind,i)=NaN;
        %Make sure fval is small (0.35 is an arbitrary cutoff)
        ind1=find(fval(:,i)>0.35);
        Dbfit(ind1,i)=NaN;
        frames=1:1:length(meanerror);
        if exist('probeoffmark')
            for m=1:size(probeoffmark,1)
                if probeoffmark(m,1)==999
                    Dbfit(1:Marksflow(probeoffmark(m,2)),:)=NaN;
                elseif probeoffmark(m,2)==999
                    Dbfit(Marksflow(probeoffmark(m,1)):end,:)=NaN;
                else
                    Dbfit(Marksflow(probeoffmark(m,1)):Marksflow(probeoffmark(m,2)),:)=NaN;
                end
            end
        end

        cutoffmeanerror(i)=10;
        goners= abs(meanerror(:,i))>cutoffmeanerror(i);
        Dbfit(goners,i)=NaN;
        clear goners
        
        %Final plot of Dbfit
        figure,plot(Dbfit(:,i),'.-','MarkerSize',20,'LineWidth',3)
        xlabel('Frame','FontSize',25)
        ylabel('BFI','FontSize',25)
        ylim([0 5e-7])
        axis tight
        tmplim=get(gca,'YLim');
        for kkkk=1:length(Marksflow)
            h=line([Marksflow(kkkk) Marksflow(kkkk)],[tmplim(1) tmplim(2)]);
            set(h,'Color',[0 0 0]);
        end
        set(gca,'FontSize',20)
        set(gcf,'PaperPositionMode','Auto')
        maxwindows(gcf)
        if fitbeta==0
            saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFit_' studyID '_S' num2str(s) '_' analysisID '_CH' num2str(i) '_fixbeta.fig'],'fig')
            saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFit_' studyID '_S' num2str(s) '_' analysisID '_CH' num2str(i) '_fixbeta.jpg'],'jpg')
        else
            saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFit_' studyID '_S' num2str(s) '_' analysisID '_CH' num2str(i) '_fitbeta.fig'],'fig')
            saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFit_' studyID '_S' num2str(s) '_' analysisID '_CH' num2str(i) '_fitbeta.jpg'],'jpg')
        end
    end

    muao=mua;
    ff=['save ' studyID '_S' num2str(s) '_' analysisID '_flow_output_fitindiv.mat timeaxis_flow mua muao musp muspo taus Dbfit Curvefit Curvefitg2 errorfit meanerror stderror corrs g1 intensitydata Marksflow Betasave numframestoavg fval exitflag usedflowdets fitbeta fitavg startcorr datalength avgnum cutoff prefix regionmarks regionlabels'];
    eval(ff);
else
    %Make sure fit converged
    ind=find(exitflagavg==0);
    Dbfitavg(ind)=NaN;
    %Make sure fval is small
    ind1=find(fvalavg>0.35);
    Dbfitavg(ind1)=NaN;

    %Plot errors in fits
    frames=1:1:length(meanerroravg);
    figure,
    [AX,h1,h2]=plotyy(frames,Dbfitavg,frames,meanerroravg);%,'b.-','MarkerSize',20,'LineWidth',3)
    set(get(AX(1),'Ylabel'),'String','BFI','FontSize',25)
    set(get(AX(2),'Ylabel'),'String','Mean Error (%)','FontSize',25)
    set(AX(2),'YTick',-50:10:50,'YLim',[-50 50],'FontSize',20)
    grid on
    set(AX(1),'YLim',[min(Dbfitavg) max(Dbfitavg)],'FontSize',20)
    set(h1,'LineWidth',3, 'LineStyle','.', 'MarkerSize', 20)
    set(h2,'LineWidth',3)
    xlabel('Frame','FontSize',25)
    title('Average of all Detectors','FontSize',30)
    set(gca,'FontSize',20)
    set(gcf,'PaperPositionMode','Auto')
    maxwindows(gcf)
    axis tight
    
    
        %If notice a correlation between erroneous Db values and mean error,
        %define cutoff for mean error
        errorquartiles = quantile(abs(nonzeros(meanerroravg)), 99);
        if cutoffmeanerror_perc
            cutoffmeanerror = errorquartiles(100-round(cutoffmeanerror_perc));
        else
            cutoffmeanerror = max(abs(nonzeros(meanerroravg)));
        end

        cutoffmeanerror=min([cutoffmeanerror 10])            

        goners=find(abs(meanerroravg)>cutoffmeanerror);
        Dbfitavgtemp  = Dbfitavg;
        Dbfitavgtemp(goners)=NaN;
        [AX,h1,h2]=plotyy(frames,Dbfitavgtemp,frames,meanerroravg);
        set(AX(2),'YTick',-50:10:50,'YLim',[-50 50],'FontSize',20)
        grid on
        set(AX(1),'FontSize',20)
        set(h1,'LineWidth',3, 'LineStyle','.', 'MarkerSize', 20)
        set(h2,'LineWidth',3)
        xlabel('Frame','FontSize',25)
        title('Average of all Detectors','FontSize',30)
        set(gca,'FontSize',20)
        set(gcf,'PaperPositionMode','Auto')
        maxwindows(gcf)
        axis tight

        Dbfitavg  = Dbfitavgtemp;
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'errorinfit_avg_' studyID '_S' num2str(s) '_' analysisID '.fig'],'fig')
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'errorinfit_avg_' studyID '_S' num2str(s) '_' analysisID '.jpg'],'jpg')


    %Final plot of Dbfit
    figure,plot(Dbfitavg,'.-','MarkerSize',20,'LineWidth',3)
    xlabel('Frame','FontSize',25)
    ylabel('BFI','FontSize',25)
    axis tight
    ylim([min(Dbfitavg) max(Dbfitavg)]);
    tmplim=get(gca,'YLim');
    for kkkk=1:length(Marksflow)
        h=line([Marksflow(kkkk) Marksflow(kkkk)],[tmplim(1) tmplim(2)]);
        set(h,'Color',[0 0 0]);
    end
    set(gca,'FontSize',20)
    maxwindows(gcf)
    set(gcf,'PaperPositionMode','Auto')
    if fitbeta==0
        
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitavg_' studyID '_S' num2str(s) '_' analysisID '_fixbeta.fig'],'fig')
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitavg_' studyID '_S' num2str(s) '_' analysisID '_fixbeta.jpg'],'jpg')
        
    else
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitavg_' studyID '_S' num2str(s) '_' analysisID '_fitbeta.fig'],'fig')
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'DbFitavg_' studyID '_S' num2str(s) '_' analysisID '_fitbeta.jpg'],'jpg')
    end
    muao=mua;
    ff=['save ' studyID '_S' num2str(s) '_' analysisID '_flow_output_fitavg.mat timeaxis_flow muao muspo taus Dbfitavg Curvefitavg Curvefitg2avg corrsavg g1avg intensityavg Marksflow numframestoavg fvalavg exitflagavg Betasaveavg usedflowdets fitbeta fitavg startcorr datalength avgnum cutoff prefix regionmarks regionlabels'];
    eval(ff);

end

if fitbeta==1
    if fitavg==0
        figure,plot(Betasave,'.','MarkerSize',10)
    elseif fitavg==1
        figure,plot(Betasaveavg,'.','MarkerSize',10)
    end
    xlabel('Frame','FontSize',25)
    ylabel('\beta','FontSize',25)
    axis tight
    ylim([0 1])
    set(gca,'FontSize',20)
    set(gcf,'PaperPositionMode','Auto')
    if fitbeta==1 & fitavg==0
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'Beta_' studyID '_S' num2str(s) '_' analysisID '.fig'],'fig')
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'Beta_' studyID '_S' num2str(s) '_' analysisID '.jpg'],'jpg')

    elseif fitbeta==1 & fitavg==1
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'Betaavg_' studyID '_S' num2str(s) '_' analysisID '.fig'],'fig')
        saveas(gcf,['..' filesep studyID filesep 'savedfigs' filesep 'Betaavg_' studyID '_S' num2str(s) '_' analysisID '.jpg'],'jpg')

    end
end
