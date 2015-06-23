close all
clear all

studyID='AUG095_032415';
analysisID='BFI_032515abs_fitavg';
s=1; % source selection, 1 or 2;
ext='';

% Remove all data when diffusion sequence occurs.  Enter start/stop marks
excludemarks = [];

%% CHECK DETECTOR AVERAGE FITS
checkingfits=1;
while checkingfits==1
    checkCH=input('Check average or individual fits? (0=avg, 1-8=individual det CH): ');
    if ~checkCH
        fname1=[ studyID '_S' num2str(s)];
        if ~isempty(analysisID)
            fname1 = [fname1 '_' analysisID];
        end
        load([ fname1 '_flow_output_fitavg.mat']);
        if fitbeta==0
            h=open(['../' studyID '/savedfigs/DbFitavg_' fname1 '_fixbeta.fig']);
        else
            h=open(['../' studyID '/savedfigs/DbFitavg_' fname1 '_fitbeta.fig']);
        end
        while checkingfits==1
            try 
                h1=get(h,'Children'); %get the axes, 3 = flow, 4= oxy
                ['Click on Line to Analyze and press a key']
                pause;

                ['Click on Point to see Correlation Curve and Fit.']

                [x,y,r]=selectpoints(gco,1,0);


                figure,semilogx(taus,corrsavg(x,:),'b-','LineWidth',2)
                hold on,semilogx(taus,Curvefitg2avg(:,x),'k--','LineWidth',3)
                ylim([0.9 1.6])
                xlim([0 1e-2])
                text(5e-4,1.3,['I=' num2str(intensityavg(x),'%6.2f') ],'FontSize',16 )
                text(5e-4,1.4,['\beta=' num2str(Betasaveavg(x),'%6.2f') ],'FontSize',16 )
                xlabel('\tau','FontSize',20)
                ylabel('g2','FontSize',20)
                set(gca,'FontSize',20)
                set(findall(gcf,'-property','FontSize'),'FontSize',20)
                set(gcf,'PaperPositionMode','Auto')
                saveas(gcf,['../' studyID '/savedfigs/fit_taupower_' fname1 '_' num2str(x) '.fig'],'fig')
                saveas(gcf,['../' studyID '/savedfigs/fit_taupower_' fname1 '_' num2str(x) '.jpg'],'jpg')

                likefit=input('Satisfied with this fit? 0=no, 1=yes  ');

                if likefit==0
                    Dbfitavg(x)=NaN;
                end
            catch
            end
            checkingfits=input('Want to look at more fits? 0=no, 1=yes  ');
        end
        close all

        if ~isempty(excludemarks)
            % Remove all data in exclude start/stop marks
            Dbfitavg(Marksflow(excludemarks(1)):Marksflow(excludemarks(2)))=NaN;
        end


        %% Final plot of Dbfit
        figure,plot(Dbfitavg,'.-','MarkerSize',20,'LineWidth',3)
        xlabel('Frame','FontSize',25)
        ylabel('BFI','FontSize',25)
        axis tight
        %ylim([0 1e-7])
        tmplim=get(gca,'YLim');
        labelshift=1;
        for kkkk=1:length(Marksflow)
            h=line([Marksflow(kkkk) Marksflow(kkkk)],[tmplim(1) tmplim(2)]);
            set(h,'Color',[0 0 0]);
        end
        for kkkk=1:length(regionlabels)
            ht=text(Marksflow(regionmarks{kkkk}(1))+labelshift, tmplim(2),regionlabels(kkkk),'Color',[0 0 0],'FontSize',16);
            set(ht,'Rotation',-90)
        end
        set(gca,'FontSize',20)
        maxwindows(gcf)
        set(gcf,'PaperPositionMode','Auto')
        if fitbeta==0
            saveas(gcf,['../' studyID '/savedfigs/DbFitavg_' fname1 '_fixbeta.fig'],'fig')
            saveas(gcf,['../' studyID '/savedfigs/DbFitavg_' fname1 '_fixbeta.jpg'],'jpg')
        else
            saveas(gcf,['../' studyID '/savedfigs/DbFitavg_' fname1 '_fitbeta.fig'],'fig')
            saveas(gcf,['../' studyID '/savedfigs/DbFitavg_' fname1 '_fitbeta.jpg'],'jpg')
        end

        ff=['save ' fname1 '_flow_output_fitavg.mat timeaxis_flow muao muspo taus Dbfitavg Curvefitavg Curvefitg2avg corrsavg g1avg intensityavg Marksflow numframestoavg fvalavg exitflagavg Betasaveavg usedflowdets fitbeta fitavg startcorr datalength avgnum cutoff prefix regionmarks regionlabels'];
        eval(ff);
    else
        %% CHECK FITS OF INDIVIDUAL PLOTS
        fname1=[ studyID '_S' num2str(s)];
        if ~isempty(analysisID)
            fname1 = [fname1 '_' analysisID];
        end
        load([ fname1 '_flow_output_fitindiv.mat']);
        if fitbeta==0
            h=open(['../' studyID '/savedfigs/DbFit_' fname1 '_CH' num2str(checkCH) '_fixbeta.fig']);
        else
            h=open(['../' studyID '/savedfigs/DbFit_' fname1 '_CH' num2str(checkCH) '_fitbeta.fig']);
        end

        xlim([Marksflow(1) Marksflow(end)])
        handles = findobj(gcf);
        
        while checkingfits==1
            try
                h1=get(h,'Children'); %get the axes, 3 = flow, 4= oxy
                ['Click on Line to Analyze and press a key']
                pause;

                ['Click on Point to see Correlation Curve and Fit.']

                [x,y,r]=selectpoints(gco,1,0);


                figure,semilogx(taus,corrs(x,:,checkCH),'b-','LineWidth',2)
                hold on,semilogx(taus,Curvefitg2(:,x,checkCH),'k--','LineWidth',3)
                ylim([0.9 1.6])
                xlim([0 1e-2])
                text(5e-4,1.3,['I=' num2str(intensitydata(x,checkCH),'%6.2f') ],'FontSize',16 )
                text(5e-4,1.4,['\beta=' num2str(Betasave(x,checkCH),'%6.2f') ],'FontSize',16 )
                xlabel('\tau','FontSize',20)
                ylabel(['g2 (CH ' num2str(checkCH) ')'],'FontSize',20)
                title(['Frame ' num2str(x)],'FontSize',20)
                set(gca,'FontSize',20)
                set(findall(gcf,'-property','FontSize'),'FontSize',20)
                set(gcf,'PaperPositionMode','Auto')
                saveas(gcf,['../' studyID '/savedfigs/fit_taupower_' fname1 '_CH' num2str(checkCH) '_' num2str(x) '.fig'],'fig')
                saveas(gcf,['../' studyID '/savedfigs/fit_taupower_' fname1 '_CH' num2str(checkCH) '_' num2str(x) '.jpg'],'jpg')

                likefit=input('Satisfied with this fit? 0=no, 1=yes  ');

                if likefit==0
                    Dbfit(x,checkCH)=NaN;
                end
            catch
            end
            
            checkingfits=input('Want to look at more fits? 0=no, 1=yes  ');
        end
        close all

        if ~isempty(excludemarks)
            % Remove all data in exclude start/stop marks
            Dbfit(Marksflow(excludemarks(1)):Marksflow(excludemarks(2)),studyID)=NaN;
        end

    end
    checkingfits = input('Check another channel? (0=no, 1=yes)');
end

if checkCH
    for checkCH = usedflowdets
        %% Final plot of Dbfit
        figure,plot(Dbfit(:,checkCH),'.-','MarkerSize',20,'LineWidth',3)
        xlabel('Frame','FontSize',25)
        ylabel(['BFI (CH ' num2str(checkCH) ')'],'FontSize',25)
        axis tight
        %ylim([0 1e-7])
        tmplim=get(gca,'YLim');
        labelshift=1;
        for kkkk=1:length(Marksflow)
            h=line([Marksflow(kkkk) Marksflow(kkkk)],[tmplim(1) tmplim(2)]);
            set(h,'Color',[0 0 0]);
        end
        for kkkk=1:length(regionlabels)
            ht=text(Marksflow(regionmarks{kkkk}(1))+labelshift, tmplim(2),regionlabels(kkkk),'Color',[0 0 0],'FontSize',16);
            set(ht,'Rotation',-90)
        end
        set(gca,'FontSize',20)
        maxwindows(gcf)
        set(gcf,'PaperPositionMode','Auto')
        xlim([Marksflow(1) Marksflow(end)])
        if fitbeta==0
            saveas(gcf,['../' studyID '/savedfigs/DbFit_' fname1 '_CH' num2str(checkCH) '_fixbeta.fig'],'fig')
            saveas(gcf,['../' studyID '/savedfigs/DbFit_' fname1 '_CH' num2str(checkCH) '_fixbeta.jpg'],'jpg')
        else
            saveas(gcf,['../' studyID '/savedfigs/DbFit_' fname1 '_CH' num2str(checkCH) '_fitbeta.fig'],'fig')
            saveas(gcf,['../' studyID '/savedfigs/DbFit_' fname1 '_CH' num2str(checkCH) '_fitbeta.jpg'],'jpg')
        end
    end

    ff=['save ' fname1 '_flow_output_fitindiv.mat timeaxis_flow muao muspo taus Dbfit Curvefit Curvefitg2 corrs g1 intensitydata Marksflow numframestoavg fval exitflag Betasave usedflowdets fitbeta fitavg startcorr datalength avgnum cutoff prefix regionmarks regionlabels'];
    eval(ff);
end
