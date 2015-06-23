%Assume 75% water
close all
clear all

studyID='AUG099_042415'; %UDPATE to match study data directory
analysisID='absolutes_042815!'; %UPDATE to match FD-DOS data directory
ext='';

waterconc=0.75;%Assume 75% water

left=[];
right=[];

fname=[studyID '_' analysisID ext '_absolute_dpfout.mat'];
load([ fname ])

left='';
right='';
leftparietal='';
rightparietal='';
calib='';
check='';
for j=1:length(grouplabels)
    tmp_left=~isempty(findstr(char(grouplabels{j}),'eft')) && ~isempty(findstr(char(grouplabels{j}),'orehead'));
    if tmp_left
        left=groups{j};
    end
    tmp_right=~isempty(findstr(char(grouplabels{j}),'ight')) && ~isempty(findstr(char(grouplabels{j}),'orehead'));
    if tmp_right
        right=groups{j};
    end
    tmp_calib=findstr(char(grouplabels{j}),'alibration');
    if ~isempty(tmp_calib)
        calib=groups{j};
    end
    tmp_check=findstr(char(grouplabels{j}),'heck');
    if ~isempty(tmp_check)
        check=groups{j};
    end
    tmp_leftparietal=~isempty(findstr(char(grouplabels{j}),'eft')) && ~isempty(findstr(char(grouplabels{j}),'ietal'));
    if tmp_leftparietal
        leftparietal=groups{j};
    end
    tmp_rightparietal=~isempty(findstr(char(grouplabels{j}),'ight')) && ~isempty(findstr(char(grouplabels{j}),'ietal'));
    if tmp_rightparietal
        rightparietal=groups{j};
    end
end

muawater=[0.0047 0.0288 0.0212];%mua h2o 688, 826nm, 786nm
%Adjust mua for water contribution
if ~isempty(left)
    for j=1:length(left)
        leftmua(j,:)=fitmua(left(j),:)-muawater(1:2)*waterconc;
    end
    leftmusp=fitmusp(left,:);
else
    leftmua=NaN(length(right),3);
    leftmusp=NaN(length(right),3);
    hbtmp_left=NaN(length(right),2);
end

if ~isempty(right)
    for j=1:length(right)
        rightmua(j,:)=fitmua(right(j),:)-muawater(1:2)*waterconc;
    end
    rightmusp=fitmusp(right,:);
else
    rightmua=NaN(length(left),3);
    rightmusp=NaN(length(left),3);
    hbtmp_right=NaN(length(left),2);
end

if ~isempty(leftparietal)
    for j=1:length(leftparietal)
        leftparietalmua(j,:)=fitmua(leftparietal(j),:)-muawater(1:2)*waterconc;
    end
    leftparietalmusp=fitmusp(leftparietal,:);
else
    leftparietalmua=NaN(length(rightparietal),3);
    leftparietalmusp=NaN(length(rightparietal),3);
    hbtmp_leftparietal=NaN(length(rightparietal),2);
end

if ~isempty(rightparietal)
    for j=1:length(rightparietal)
        rightparietalmua(j,:)=fitmua(rightparietal(j),:)-muawater(1:2)*waterconc;
    end
    rightparietalmusp=fitmusp(rightparietal,:);
else
    rightparietalmua=NaN(length(leftparietal),3);
    rightparietalmusp=NaN(length(leftparietal),3);
    hbtmp_rightparietal=NaN(length(leftparietal),2);
end

calibmua=nanmean(fitmua(calib,:),1);
checkmua=nanmean(fitmua(check,:),1);
calibmusp=nanmean(fitmusp(calib,:),1);
checkmusp=nanmean(fitmusp(check,:),1);

ext1=[276 2051.96; 974 693.04]*2.303/1e6;%W. B. Gratzer, converted to uM; from Prahl's website
ext786=[740	957.36]*2.303/1e6;%W. B. Gratzer, converted to uM; from Prahl's website [HbO2(786) Hb(786)]

for i=1:size(rightmua)
    
    hbtmp_right(i,:)=inv(ext1)*rightmua(i,1:2).';
   
    Hb_right(i)=hbtmp_right(i,2);
    HbO2_right(i)=hbtmp_right(i,1);
    
    THC_right(i)=sum(hbtmp_right(i,:));
    StO2_right(i)=hbtmp_right(i,1)./THC_right(i)*100;
    
    %Calculate mua at 785 from extinction coef's and HbO2/Hb concentrations
    rightmua(i,3)=[HbO2_right(i) Hb_right(i)]*ext786.'+muawater(3)*waterconc;
    %Calculate musp at 785 by averaging 685 and 830 (crude, I know...)
    rightmusp(i,3)=(rightmusp(i,1)+rightmusp(i,2))/2;
    DPF_right(i,:)=3.*rightmusp(i,:).*2.5./2./(2.5*sqrt(3*rightmua(i,:).*rightmusp(i,:))+1);

end

for i=1:size(rightparietalmua)
    
    hbtmp_rightparietal(i,:)=inv(ext1)*rightparietalmua(i,1:2).';
   
    Hb_rightparietal(i)=hbtmp_rightparietal(i,2);
    HbO2_rightparietal(i)=hbtmp_rightparietal(i,1);
    
    THC_rightparietal(i)=sum(hbtmp_rightparietal(i,:));
    StO2_rightparietal(i)=hbtmp_rightparietal(i,1)./THC_rightparietal(i)*100;
    
    %Calculate mua at 785 from extinction coef's and HbO2/Hb concentrations
    rightparietalmua(i,3)=[HbO2_rightparietal(i) Hb_rightparietal(i)]*ext786.'+muawater(3)*waterconc;
    %Calculate musp at 785 by averaging 685 and 830 (crude, I know...)
    rightparietalmusp(i,3)=(rightparietalmusp(i,1)+rightparietalmusp(i,2))/2;
    DPF_rightparietal(i,:)=3.*rightparietalmusp(i,:).*2.5./2./(2.5*sqrt(3*rightparietalmua(i,:).*rightparietalmusp(i,:))+1);

end

for i=1:size(leftmua)
    
    hbtmp_left(i,:)=inv(ext1)*leftmua(i,1:2).';
    
    Hb_left(i)=hbtmp_left(i,2);
    HbO2_left(i)=hbtmp_left(i,1);
    
    THC_left(i)=sum(hbtmp_left(i,:));
    StO2_left(i)=hbtmp_left(i,1)./THC_left(i)*100;
    
    %Calculate mua at 785 from extinction coef's and HbO2/Hb concentrations
    leftmua(i,3)=[HbO2_left(i) Hb_left(i)]*ext786.'+muawater(3)*waterconc;
    %Calculate musp at 785 by averaging 685 and 830 (crude, I know...)
    leftmusp(i,3)=(leftmusp(i,1)+leftmusp(i,2))/2;
    DPF_left(i,:)=3.*leftmusp(i,:).*2.5./2./(2.5*sqrt(3*leftmua(i,:).*leftmusp(i,:))+1);
    
end

for i=1:size(leftparietalmua)
    
    hbtmp_leftparietal(i,:)=inv(ext1)*leftparietalmua(i,1:2).';
    
    Hb_leftparietal(i)=hbtmp_leftparietal(i,2);
    HbO2_leftparietal(i)=hbtmp_leftparietal(i,1);
    
    THC_leftparietal(i)=sum(hbtmp_leftparietal(i,:));
    StO2_leftparietal(i)=hbtmp_leftparietal(i,1)./THC_leftparietal(i)*100;
    
    %Calculate mua at 785 from extinction coef's and HbO2/Hb concentrations
    leftparietalmua(i,3)=[HbO2_leftparietal(i) Hb_leftparietal(i)]*ext786.'+muawater(3)*waterconc;
    %Calculate musp at 785 by averaging 685 and 830 (crude, I know...)
    leftparietalmusp(i,3)=(leftparietalmusp(i,1)+leftparietalmusp(i,2))/2;
    DPF_leftparietal(i,:)=3.*leftparietalmusp(i,:).*2.5./2./(2.5*sqrt(3*leftparietalmua(i,:).*leftparietalmusp(i,:))+1);
    
end

DPF_left=nanmean(DPF_left,1);
DPF_right=nanmean(DPF_right,1);
% DPF_leftparietal=nanmean(DPF_leftparietal,1);
% DPF_rightparietal=nanmean(DPF_rightparietal,1);

figure,subplot(1,2,1)
plot(StO2_left,'.','MarkerSize',40,'LineWidth',3)
hold on,plot(StO2_right,'o','MarkerSize',15,'LineWidth',3)
legend('Left Forehead','Right Forehead')
ylabel('StO2(%)')
xlim([0.5 length(StO2_right)+0.5])
set(gca,'XTick',1:length(THC_right))
%set(gca,'XTick',t)
xlabel('Repetition No.')
grid on
ylim([10 80])

subplot(1,2,2)
plot(THC_left,'.','MarkerSize',40,'Color',[1 0.6 0],'LineWidth',3)
hold on,plot(Hb_left,'b.','MarkerSize',40,'LineWidth',3)
hold on,plot(HbO2_left,'r.','MarkerSize',40,'LineWidth',3)

hold on,plot(THC_right,'o','Color',[1 0.6 0],'MarkerSize',15,'LineWidth',3)
hold on,plot(Hb_right,'bo','MarkerSize',15,'LineWidth',3)
hold on,plot(HbO2_right,'ro','MarkerSize',15,'LineWidth',3)

xlim([0.5 length(THC_right)+0.5])
ylim([0 150])
set(gca,'XTick',1:length(THC_right))
ylabel('\mu M ')
legend('THC','Hb','HbO2')
xlabel('Repetition No.')
grid on

maxwindows(gcf)
set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(gcf,'PaperPositionMode','Auto')
saveas(gcf,['../' studyID '/savedfigs/Absoluteoxygenation_' studyID '_' analysisID ext '.fig'],'fig')
saveas(gcf,['../' studyID '/savedfigs/Absoluteoxygenation_' studyID '_' analysisID ext '.jpg'],'jpg')


if ~isempty(strfind(analysisID,'!')) && (exist('StO2_leftparietal','var') || exist('StO2_rightparietal','var'))
    figure,subplot(1,2,1)
    plot(StO2_leftparietal,'*','MarkerSize',15,'LineWidth',3)
    hold on,plot(StO2_rightparietal,'+','MarkerSize',15,'LineWidth',3)
    legend('Left Parietal','Right Parietal')
    ylabel('StO2(%)')
    xlim([0.5 length(StO2_right)+0.5])
    set(gca,'XTick',1:length(THC_rightparietal))
    %set(gca,'XTick',t)
    xlabel('Repetition No.')
    grid on
    ylim([10 80])

    subplot(1,2,2)
    plot(THC_leftparietal,'*','MarkerSize',15,'Color',[1 0.6 0],'LineWidth',3)
    hold on,plot(Hb_leftparietal,'b*','MarkerSize',15,'LineWidth',3)
    hold on,plot(HbO2_leftparietal,'r*','MarkerSize',15,'LineWidth',3)

    hold on,plot(THC_rightparietal,'+','Color',[1 0.6 0],'MarkerSize',15,'LineWidth',3)
    hold on,plot(Hb_rightparietal,'b+','MarkerSize',15,'LineWidth',3)
    hold on,plot(HbO2_rightparietal,'r+','MarkerSize',15,'LineWidth',3)

    xlim([0.5 length(THC_right)+0.5])
    ylim([0 150])
    set(gca,'XTick',1:length(THC_right))
    ylabel('\mu M ')
    legend('THC','Hb','HbO2')
    xlabel('Repetition No.')
    grid on
end

maxwindows(gcf)
set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(gcf,'PaperPositionMode','Auto')
saveas(gcf,['../' studyID '/savedfigs/Absoluteoxygenation_parietal_' studyID '_' analysisID ext '.fig'],'fig')
saveas(gcf,['../' studyID '/savedfigs/Absoluteoxygenation_parietal_' studyID '_' analysisID ext '.jpg'],'jpg')


maxwindows(gcf)
set(findall(gcf,'-property','FontSize'),'FontSize',30)
set(gcf,'PaperPositionMode','Auto')
saveas(gcf,['../' studyID '/savedfigs/Absoluteoxygenation_' studyID '_' analysisID ext '.fig'],'fig')
saveas(gcf,['../' studyID '/savedfigs/Absoluteoxygenation_' studyID '_' analysisID ext '.jpg'],'jpg')


if ~(exist('StO2_leftparietal','var') || exist('StO2_rightparietal','var')) 
    tt=['save ' studyID '_' analysisID ext '_baselines.mat leftmua rightmua leftmusp rightmusp HbO2_left HbO2_right Hb_left Hb_right THC_left THC_right StO2_left StO2_right DPF_left DPF_right'];
    ttt=cat(2,nanmean(StO2_right),nanstd(StO2_right),nanmean(StO2_left),nanstd(StO2_left),nanmean(THC_right),nanstd(THC_right),nanmean(THC_left),nanstd(THC_left));
else
    tt=['save ' studyID '_' analysisID ext '_baselines.mat leftmua rightmua leftparietalmua rightparietalmua leftmusp rightmusp HbO2_left HbO2_right Hb_left Hb_right THC_left THC_right StO2_left StO2_right DPF_left DPF_right leftparietalmusp rightparietalmusp HbO2_leftparietal HbO2_rightparietal Hb_leftparietal Hb_rightparietal THC_leftparietal THC_rightparietal StO2_leftparietal StO2_rightparietal DPF_leftparietal DPF_rightparietal'];
    ttt=cat(2,nanmean(StO2_right),nanstd(StO2_right),nanmean(StO2_left),nanstd(StO2_left),nanmean(StO2_rightparietal),nanstd(StO2_rightparietal),nanmean(StO2_leftparietal),nanstd(StO2_leftparietal),nanmean(THC_right),nanstd(THC_right),nanmean(THC_left),nanstd(THC_left),nanmean(THC_rightparietal),nanstd(THC_rightparietal),nanmean(THC_leftparietal),nanstd(THC_leftparietal));ttt=cat(2,nanmean(StO2_right),nanstd(StO2_right),nanmean(StO2_left),nanstd(StO2_left),nanmean(StO2_rightparietal),nanstd(StO2_rightparietal),nanmean(StO2_leftparietal),nanstd(StO2_leftparietal),nanmean(THC_right),nanstd(THC_right),nanmean(THC_left),nanstd(THC_left),nanmean(THC_rightparietal),nanstd(THC_rightparietal),nanmean(THC_leftparietal),nanstd(THC_leftparietal));
end
eval(tt);
cat(1,ttt,nan(size(ttt)))

muaforxl=nanmean(cat(1,leftmua(:,1:2),rightmua(:,1:2)),1);
muspforxl=nanmean(cat(1,leftmusp(:,1:2),rightmusp(:,1:2)),1);
ttt=cat(2,muaforxl,muspforxl);
cat(1,ttt,nan(size(ttt)))

allHbO2=cat(2,HbO2_right,HbO2_left);
allHb=cat(2,Hb_right,Hb_left);

ttt=cat(2,nanmean(allHb),nanstd(allHb),nanmean(allHbO2),nanstd(allHbO2));
cat(1,ttt,nan(size(ttt)))
