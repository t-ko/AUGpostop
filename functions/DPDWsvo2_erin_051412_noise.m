%Updated 7/6/09 to subtract offsets from I and Q data!  We are very close
%to noise floor, so it is important to do this.
%Updated 7/6/09 to include calculation of mua changes and to work for both
%ISS and Joel's Instruments

function DPDWsvo2_erin_051412_1(noiselevel,ext)

patID='erin';
patdate='051412';
probeoffmark=[];
baselinemarks=[];
age=27;
%If you want something added to your file name
%Extension you want added to filename, if nothing, just do ext='';
%ext=extension;
filenum=1;
%syncmark=1;
appliedcalibration=1;

%infant DPF
%DPF830=4.35;
%DPF690=4.2;
%DPF790=4.4;
DPF830=4.67 + 0.62*age^0.819;
DPF690=5.38 + 0.049*age^0.877;
DPF790=4.99 + 0.067*age^0.814;

%DPFvalues from Duncan et al. 1995 paper in Phys. Med. Biol.  Optical
%pathlength measurements on adult head, calf...


usedISS=1;%0=Joel's Instrument, 1=ISS
usednewsettings=0;%Dennis made new settings file for us to turn off PMTs while DCS laser is on.  
useddet=2;%If using ISS instrument, DetA=1, DetB=2
threelambda=0;%Sometimes we have trouble with 685 nm laser, if this is the case, just use 785 and 830 for DPF calcs

%DPF Values From Erin

DPF=[DPF690 DPF830];
framesperfile=1;

if usedISS
    %ISS Analysis
    sourcesused=[13 12]; 
    sourcesnotused=[1 2];
    %r=[1.44 1.94 2.41 2.90; 1.51 1.99 2.49 2.98];
    r=[2.5; 2.5];
    sourceindex=[1;2];

    lambdas=[688 826];

    %LOAD DATA
    fdir=['/home/jlynch/Data/CHOP/PH/svo2test/svo2test_' patdate '/' patID '/'];
    % Must sort the rest
    files=dir([ fdir '*.txt']);
    for numfile=1:size(files,1)
        names(numfile,:)=files(numfile).name(12:end-4);
    end
    names=sortrows(names);

    data=[];
    %for i=1:size(names,1)
    i=filenum;
        %Find time each file was saved "HH:MM:SS"
        timetmp{1}=[ num2str(names(i,2:3)) ':' num2str(names(i,4:5)) ':' num2str(names(i,6:7)) ];
        fname=['/home/jlynch/Data/CHOP/PH/svo2test/svo2test_' patdate '/' patID '/Data_' patdate names(i,:) '.txt'];
        fid = fopen([ fname ], 'r');
        if usednewsettings
            [tmpdata, count]=fscanf(fid,'%c %s',1010);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%g %g',[79 inf]);
 elseif appliedcalibration
            [tmpdata, count]=fscanf(fid,'%c %s',657);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%c %s',96);
            phacalibcoeffs=str2num(tmpdata);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%c %s',2);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%c %s',64);
            ACDCcalibcoeffs=str2num(tmpdata);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%c %s',390);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%g %g',[103 inf]);
    else
            [tmpdata, count]=fscanf(fid,'%c %s',1210);
            clear tmpdata count
            [tmpdata, count]=fscanf(fid,'%g %g',[103 inf]);
    end
        fclose(fid);
        tmpdata=tmpdata.';
        %Calculate how many frames were recorded in the given file
        datalength(1)=size(tmpdata,1)/framesperfile;
        
        data=cat(1,data,tmpdata);

    %end

    % The columns of data:
    % Column 1: Time (in seconds) *NOT RELIABLE! DON'T USE FOR TIME AXIS!*
    % Column 2: Data group
    % Column 3: Step (has to do with external GPIB devices)
    % Column 4: Mark
    % Column 5: Flag (turn into binary for message about data)
    % Columns 6-21: AC data for Detector A, Sources 1-16
    % Columns 22-37: DC data for Detector A, Sources 1-16
    % Columns 38-53: Phase data for zDetector A, Sources 1-16
    % Columns 54-69: AC data for Detector B, Sources 1-16
    % Columns 70-85: DC data for Detector B, Sources 1-16
    % Columns 86-101: Phase data for Detector B, Sources 1-16
    % Columns 102: Analog output data 1
    % Column 103: Digital output data

    if usednewsettings
        updatetime=1/6.3776;
    else
        updatetime=1/6.25; % Update time, set to 6.25 Hz.
    end
    Marks=find(data(:,4));
    Marks=floor(Marks/framesperfile);

    % Plot noncalibrated AC, DC and phase data for Detectors A & B
    for j=1:length(sourcesused)
        if useddet==1
            ind=5;
        elseif useddet==2
            ind=53;
        end
          if appliedcalibration
        ACtmp(:,j)=data(:,ind+sourcesused(j))/ACDCcalibcoeffs(sourcesused(j)) - noiselevel*data(:,ind+sourcesused(j))/ACDCcalibcoeffs(sourcesused(j));
        DCtmp(:,j)=data(:,ind+16+sourcesused(j))/ACDCcalibcoeffs(16+sourcesused(j))-noiselevel*data(:,ind+16+sourcesused(j))/ACDCcalibcoeffs(16+sourcesused(j));
        phasetmp(:,j)=(data(:,ind+16*2+sourcesused(j))+phacalibcoeffs(32+sourcesused(j)))*pi/180-noiselevel*(data(:,ind+16*2+sourcesused(j))+phacalibcoeffs(32+sourcesused(j)))*pi/180;
    else
        ACtmp(:,j)=data(:,ind+sourcesused(j))-noiselevel*data(:,ind+sourcesused(j));
        DCtmp(:,j)=data(:,ind+16+sourcesused(j))-noiselevel*data(:,ind+16+sourcesused(j));
        phasetmp(:,j)=data(:,ind+16*2+sourcesused(j))-noiselevel*data(:,ind+16*2+sourcesused(j));
          end
        %Convert phase to radians
        phasetmp(:,j)=phasetmp(:,j)*pi/180;
        phasetmp(:,j)=unwrap(phasetmp(:,j));
        %Convert back to degrees
        phasetmp(:,j)=phasetmp(:,j)*180/pi;
        
        %Get data from sources 2,3, and 5
        ACtmp_offsource(:,j)=data(:,ind+sourcesnotused(j));
        DCtmp_offsource(:,j)=data(:,ind+16+sourcesnotused(j));
        phasetmp_offsource(:,j)=data(:,ind+16*2+sourcesnotused(j));
        %Convert phase to radians
        phasetmp_offsource(:,j)=phasetmp_offsource(:,j)*pi/180;
        phasetmp_offsource(:,j)=unwrap(phasetmp_offsource(:,j));
        %Convert back to degrees
        phasetmp_offsource(:,j)=phasetmp_offsource(:,j)*180/pi;
    end
    
    
  

    %Find the mean and standard deviation of each frame of data
    for i=1:length(data(:,2))/framesperfile
        AC(i,:)=mean(ACtmp((i-1)*framesperfile+1:framesperfile*i,:),1);
        DC(i,:)=mean(DCtmp((i-1)*framesperfile+1:framesperfile*i,:),1);
        phase(i,:)=mean(phasetmp((i-1)*framesperfile+1:framesperfile*i,:),1);
        ACstd(i,:)=std(ACtmp((i-1)*framesperfile+1:framesperfile*i,:),0,1);
        DCstd(i,:)=std(DCtmp((i-1)*framesperfile+1:framesperfile*i,:),0,1);
        phasestd(i,:)=std(phasetmp((i-1)*framesperfile+1:framesperfile*i,:),0,1);
        
        AC_offsource(i,:)=mean(ACtmp_offsource((i-1)*framesperfile+1:framesperfile*i,:),1);
        DC_offsource(i,:)=mean(DCtmp_offsource((i-1)*framesperfile+1:framesperfile*i,:),1);
        phase_offsource(i,:)=mean(phasetmp_offsource((i-1)*framesperfile+1:framesperfile*i,:),1);
        ACstd_offsource(i,:)=std(ACtmp_offsource((i-1)*framesperfile+1:framesperfile*i,:),0,1);
        DCstd_offsource(i,:)=std(DCtmp_offsource((i-1)*framesperfile+1:framesperfile*i,:),0,1);
        phasestd_offsource(i,:)=std(phasetmp_offsource((i-1)*framesperfile+1:framesperfile*i,:),0,1);
        
    end

    %GET TIME AXIS!!! %See lab notebook 4/30/10
    ISStime4frame=datenum(timetmp,'HH:MM:SS')-floor(datenum(timetmp,'HH:MM:SS'));%In arbitrary units--1a.u.=24hrs, counting from 1/1/2000.
    ISStime4frame(find(datalength==0))=[];
    datalength(find(datalength==0))=[];
    %Subtract off whole days since year 2000 and you are left with a fraction
    %of a day that has elapsed since the start of this data
    
    %Zero time axis so that t=0 is first frame
    ISStime4frame=(ISStime4frame-ISStime4frame(1))*1.44e3;%Convert time in a.u. to minutes
    
    %Now, we only have time points for when each file was saved.  However,
    %each file has approximately 5 frames.  So we will interpolate the time
    %points for the frames within each file.
    
    %Get frame number that each ISS time corresponds to:
    for i=1:length(ISStime4frame)
        frame(i+1)=sum(datalength(1:i));
    end
    %Since first time point corresponds to the first file, we need to add a fake file in the
    %beginning in order for interpolation to catch frames 1-5
    frame(1)=0;
    ISStime4frame(2:end+1)=ISStime4frame(1:end);
    %ISStime4frame(1)=ISStime4frame(2)-(ISStime4frame(3)-ISStime4frame(2));
    
    %Interpolate times for the rest of the frames
    %First find index of jumps in time (we dont want to interpolate in the
    %jumps)
    ind(1)=1;
    ind(2)=length(ISStime4frame);
    for j=1:size(ind,2)-1
        %Define chunk of continuous data without jump in time
        chunk{j}=ind(j):1:ind(j+1);
        if j>1
            ISStime4frame(ind(j))=ISStime4frame(ind(j)+1)-(ISStime4frame(ind(j)+2)-ISStime4frame(ind(j)+1));
        end
        %Interpolate time data within this chunk
        ISStime(min(frame(chunk{j}))+1:1:max(frame(chunk{j}))+1)=interp1(frame(chunk{j}),ISStime4frame(chunk{j}),min(frame(chunk{j})):1:max(frame(chunk{j})));
    end
    ISStimetmp=ISStime(2:end);%Ignore first point, just used it to interpolate 
    clear ind ISStime
    %ISStimetmp=ISStimetmp-ISStimetmp(Marks(syncmark));
    ISStime=ISStimetmp;%Now we have timepoints for each datapoint!!!
    

    %load the probe map, sdlist etc.
    plott=0;
    %flatbabyprobe;%Loads s-d seps in cm
sdlist(1,1)=2.5;
    
    
    %Every time probe moves, start new baselines
    if size(baselinemarks)==0
             baselineoxy=1:200; %baseline frames
      
            if threelambda
                [hbseries,hbo2series,muaseries]=threewavelengthdpf(squeeze(AC(:,1)),squeeze(AC(:,2)),squeeze(AC(:,3)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,2))),mean(squeeze(AC(baselineoxy,3))),...
                    lambdas,sdlist(1,1),DPF);
            else
                [hbseries,hbo2series,muaseries]=twowavelengthdpf(squeeze(AC(:,1)),squeeze(AC(:,2)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,2))),...
                    lambdas(1:2),sdlist(1,1),DPF(1:2));
            end
      
        clear baselineoxy
    else
    for k=1:size(baselinemarks,1)
        baselineoxy=Marks(baselinemarks(k,1)):Marks(baselinemarks(k,2)); %baseline frames
        if k==1
            if threelambda
                [hbseries,hbo2series,muaseries]=threewavelengthdpf(squeeze(AC(:,1)),squeeze(AC(:,2)),squeeze(AC(:,3)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,2))),mean(squeeze(AC(baselineoxy,3))),...
                    lambdas,sdlist(1,1),DPF);
            else
                [hbseries,hbo2series,muaseries]=twowavelengthdpf(squeeze(AC(:,1)),squeeze(AC(:,2)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,2))),...
                    lambdas(1:2),sdlist(1,1),DPF(1:2));
            end
        else
            if threelambda
                [hbseries(baselineoxy(1):length(AC)),hbo2series(baselineoxy(1):length(AC)),muaseries(:,baselineoxy(1):length(AC))]=threewavelengthdpf(squeeze(AC(baselineoxy(1):end,1)),squeeze(AC(baselineoxy(1):end,2)),squeeze(AC(baselineoxy(1):end,3)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,2))),mean(squeeze(AC(baselineoxy,3))),...
                    lambdas,sdlist(1,1),DPF);
            else
                [hbseries(baselineoxy(1):length(AC)),hbo2series(baselineoxy(1):length(AC)),muaseries(:,baselineoxy(1):length(AC))]=twowavelengthdpf(squeeze(AC(baselineoxy(1):end,1)),squeeze(AC(baselineoxy(1):end,3)),...
                    mean(squeeze(AC(baselineoxy,1))),mean(squeeze(AC(baselineoxy,3))),...
                    lambdas(1:2:3),sdlist(1,1),DPF(1:2:3));
            end

        end
        clear baselineoxy
    end
    end
    if exist([ patID '_' patdate '_baselines.mat'],'file')
        tt=['load ' patID '_' patdate '_baselines.mat'];
        eval(tt);
        
        ext1=[276 2051.96; 974 693.04]*2.303/1e6;%W. B. Gratzer, converted to uM; from Prahl's website
        
        for j=2:size(headmua,1)+1
            hbtmp_head(j,:)=inv(ext1)*headmua(j-1,1:2).';
            
            Hb_head(j)=hbtmp_head(j,2);
            HbO2_head(j)=hbtmp_head(j,1);
            
            THC_head(j)=sum(hbtmp_head(j,:));
            StO2_head(j)=hbtmp_head(j,1)./THC_head(j)*100;
            
        end
        %Calculate mean hemoglobins by averaging both sides of head and all
        %repetitions
        Hbmean=nanmean(Hb_head);
        Hbstd=nanstd(Hb_head);
        HbO2mean=nanmean(HbO2_head);
        HbO2std=nanstd(HbO2_head);
        %Use hemoglobin concentrations to calculate mean mua at 785
        ext2=[2051.96 276; 693.04 974; 921.8 748 ]*2.303/1e6;%W. B. Gratzer, converted to uM; from Prahl's website [830(Hb,HbO2) 690(Hb,HbO2) 788(Hb,HbO2)]
        muamean=ext2*[Hbmean;HbO2mean];
        %Adjust deltamua, dHb, dHbO2 for baseline values
        for i=1:length(lambdas)
            muaseries(i,:)=muamean(i)+muaseries(i,:);
        end
        hbseries=Hbmean+hbseries.*1000;
        hbo2series=HbO2mean+hbo2series.*1000;
        %Calculate mean musp
        muspmean=nanmean(headmusp,1);
        muspstd=nanstd(headmusp,0,1);
    else
        hbseries=hbseries.*1000;
        hbo2series=hbo2series.*1000;
    end
   
    
    %If probe was off head, set values to NaN
    for m=1:size(probeoffmark,1)
        if probeoffmark(m,1)==999
            hbseries(1:Marks(probeoffmark(m,2)))=NaN;
            hbo2series(1:Marks(probeoffmark(m,2)))=NaN;
            muaseries(:,1:Marks(probeoffmark(m,2)))=NaN;
        
        elseif probeoffmark(m,2)==999
            hbseries(Marks(probeoffmark(m,1)):end)=NaN;
            hbo2series(Marks(probeoffmark(m,1)):end)=NaN;
            muaseries(:,Marks(probeoffmark(m,1)):end)=NaN;
            
        else
            hbseries(Marks(probeoffmark(m,1)):Marks(probeoffmark(m,2)))=NaN;
            hbo2series(Marks(probeoffmark(m,1)):Marks(probeoffmark(m,2)))=NaN;
            muaseries(:,Marks(probeoffmark(m,1)):Marks(probeoffmark(m,2)))=NaN;
        end
    end
    
   
    tt=['save ' patID '_' patdate ext '_svo2_dpfout.mat ISStime AC ACstd Marks hbo2series hbseries muaseries useddet DPF usedISS threelambda probeoffmark baselinemarks'];
    eval(tt);
end