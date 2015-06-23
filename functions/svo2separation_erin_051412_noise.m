function [fwhb,fwhbo2,meansvo2,stdsvo2]=svo2separation_erin_051412_noise(ext)



patID='erin';
patdate='051412';

hr=70;
bpm=11;

percwidth=.01;
filtersize=10;

tt=['load ' patID '_' patdate ext '_svo2_dpfout.mat'];
eval(tt);
frametime=.1664; %seconds
hearttime=1/(hr/60); %seconds
breathtime=1/(bpm/60); % seconds



framerate=1/frametime;
heartrate=1/hearttime;
breathrate=1/breathtime;

timelength=size(ISStime,2);

timeaxis=[1:timelength]*frametime; % time axis in seconds



start=500;
timeaxis=timeaxis(start:end);
hbo2series=hbo2series(start:end);
hbseries=hbseries(start:end);

baseline=1; %looking at only pre-cuff data
if baseline
    baselinename='_baseline';
else
    baselinename='';
end





Y = fft(hbseries);
Y(1)=[];
n=length(Y);
power = abs(Y(1:floor(n/2))).^2;
nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist;
val1=(bpm/60)*frametime-percwidth*(bpm/60)*frametime;
val2=(bpm/60)*frametime+percwidth*(bpm/60)*frametime;
[mindiff1,f1]=(min(abs(freq-val1)));
[mindiff2,f2]=(min(abs(freq-val2)));

x=freq(f1:f2);
y=power(f1:f2);

fwhb=fwhm(x,y);


Y = fft(hbo2series);
Y(1)=[];
n=length(Y);
power = abs(Y(1:floor(n/2))).^2;
nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist;
val1=(bpm/60)*frametime-percwidth*(bpm/60)*frametime;
val2=(bpm/60)*frametime+percwidth*(bpm/60)*frametime;
[mindiff1,f1]=(min(abs(freq-val1)));
[mindiff2,f2]=(min(abs(freq-val2)));


x=freq(f1:f2);
y=power(f1:f2);

fwhbo2=fwhm(x,y);





heartrateframe=heartrate*frametime*2;
breathrateframe=breathrate*frametime*2;


[bbr,abr]=butter(1, [(breathrateframe-percwidth*breathrateframe), (breathrateframe+percwidth*breathrateframe)]);
filtbr_hbo2=filter(bbr,abr,hbo2series);
filtbr_hb=filter(bbr,abr,hbseries);

filtbr_hbo2=strokefilter(filtbr_hbo2,filtersize);
filtbr_hb=strokefilter(filtbr_hb, filtersize);


minpeakdistance_br=round((breathtime/3)/frametime);

[pks_hbo2_br,locs_pks_hbo2_br]=findpeaks(filtbr_hbo2, 'minpeakdistance',minpeakdistance_br);

[trs_hbo2_br,locs_trs_hbo2_br]=findpeaks(filtbr_hbo2*-1, 'minpeakdistance',minpeakdistance_br);


[pks_hb_br,locs_pks_hb_br]=findpeaks(filtbr_hb, 'minpeakdistance',minpeakdistance_br);

[trs_hb_br,locs_trs_hb_br]=findpeaks(filtbr_hb*-1, 'minpeakdistance',minpeakdistance_br);


if length(locs_pks_hbo2_br) < length(locs_trs_hbo2_br)
    endpeak_hbo2_br=length(locs_pks_hbo2_br)-1;
else
    endpeak_hbo2_br=length(locs_trs_hbo2_br)-1;
end

if length(locs_pks_hb_br) < length(locs_trs_hb_br)
    endpeak_hb_br=length(locs_pks_hbo2_br)-1;
else
    endpeak_hb_br=length(locs_trs_hbo2_br)-1;
end


if locs_pks_hbo2_br(1)<locs_trs_hbo2_br(1)
    for i=1:endpeak_hbo2_br-1
        hbo2maxtmp=(filtbr_hbo2(locs_pks_hbo2_br(i))+filtbr_hbo2(locs_pks_hbo2_br(i+1)))/2;
        hbo2max(i)=(filtbr_hbo2(locs_pks_hbo2_br(i))+filtbr_hbo2(locs_pks_hbo2_br(i+1)))/2;
        hbo2mintmp=filtbr_hbo2(locs_trs_hbo2_br(i));
        hbo2min(i)=filtbr_hbo2(locs_trs_hbo2_br(i));
        hbmaxtmp=(filtbr_hb(locs_pks_hb_br(i))+filtbr_hb(locs_pks_hb_br(i+1)))/2;
        hbmax(i)=hbmaxtmp;
        hbmintmp=filtbr_hb(locs_trs_hb_br(i));
        hbmin(i)=hbmintmp;
        deltahbo2tmp=hbo2maxtmp-hbo2mintmp;
        deltahbo2(i)=deltahbo2tmp;
        deltahbtmp=hbmaxtmp-hbmintmp;
        deltahb(i)=deltahbtmp;
        svo2(i)=deltahbo2tmp/(deltahbo2tmp+deltahbtmp);
        clear hbo2maxtmp;
        clear hbo2mintmp;
        clear hbmaxtmp;
        clear hbmintmp;
        clear deltahbo2tmp;
        clear deltahbtmp;
    end
else
    for i=1:endpeak_hbo2_br-1
        hbo2maxtmp=filtbr_hbo2(locs_pks_hbo2_br(i));
        hbo2mintmp=(filtbr_hbo2(locs_trs_hbo2_br(i))+filtbr_hbo2(locs_trs_hbo2_br(i+1)))/2;
        hbmaxtmp=filtbr_hb(locs_pks_hb_br(i));
        hbmintmp=(filtbr_hb(locs_trs_hb_br(i))+filtbr_hb(locs_trs_hb_br(i+1)))/2;
        deltahbo2tmp=hbo2maxtmp-hbo2mintmp;
        deltahbtmp=hbmaxtmp-hbmintmp;
        hbo2max(i)=hbo2maxtmp;
        hbo2min(i)=hbo2mintmp;
        hbmax(i)=hbmaxtmp;
        hbmin(i)=hbmintmp;
        deltahbo2(i)=deltahbo2tmp;
        deltahb(i)=deltahbtmp;
        svo2(i)=deltahbo2tmp/(deltahbo2tmp+deltahbtmp);
        clear hbo2maxtmp;
        clear hbo2mintmp;
        clear hbmaxtmp;
        clear hbmintmp;
        clear deltahbo2tmp;
        clear deltahbtmp;
    end
end

meansvo2=nanmean(svo2);
stdsvo2=std(svo2);


ff=['save ' patID '_' patdate ext '_svo2_separation_output.mat svo2 heartrate breathrate framerate'];
eval(ff);
