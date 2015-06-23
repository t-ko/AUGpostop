function [CR,linFit] = BlandAltmanSVC(var1, var2, flag)
 
    %%%Plots a Bland-Altman Plot
    %%%INPUTS:
    %%% var1 and var2 - vectors of the measurements
    %%%flag - how much you want to plot
        %%% 0 = no plot
        %%% 1 = just the data
        %%% 2 = data and the difference and CR lines
        %%% 3 = above and a linear fit
    %%%
    %%%OUTPUTS:
    %%% means = the means of the data
    %%% diffs = the raw differences
    %%% meanDiff = the mean difference
    %%% CR = the 2SD confidence limits
    %%% linfit = the paramters for the linear fit
    
    
    if (nargin<1)
       %%%Use test data
       var1=[512,430,520,428,500,600,364,380,658,445,432,626,260,477,259,350,451];%,...
       var2=[525,415,508,444,500,625,460,390,642,432,420,605,227,467,268,370,443];
       flag = 3;
    end
    
    if nargin==2
        flag = 0;
    end
    mdl = fitlm(var1,var2); %%%work out the linear fit coefficients
    linFit = [double(mdl.Coefficients(2,1)) double(mdl.Coefficients(1,1))];
    
    means = mean([var1;var2]);
    ratios = var2./var1;
    
    meanDiff = mean(ratios);
    sdDiff = std(ratios);
    CR = [1.96 * sdDiff + linFit(2), - 1.96 * sdDiff + linFit(2)]; %%95% confidence range
    p = double(mdl.Coefficients(2,end));
    R = mdl.Rsquared.Adjusted;
    %%%plot results unless flag is 0
    if flag ~= 0
        plot(var1,var2,'LineStyle','o','Color','red')
        hold on
        if flag > 1
            plot(var1, var1.*linFit(1)+CR(1),'LineStyle','--','Color','black'); %%%plot the upper CR
            plot(var1, var1.*linFit(1)+CR(2),'LineStyle','--','Color','black'); %%%plot the lower CR
        end
        if flag > 2
            plot(var1, var1.*linFit(1)+linFit(2),'LineStyle','--','Color','red'); %%%plot the linear fit            
        end        
    end
    axis tight
    text(mean(get(gca,'XLim')), mean(get(gca,'YLim')),{['y = ' num2str(linFit(1)) ' * x + ' num2str(linFit(2))];['CR: ' num2str(CR(1)) ', ' num2str(CR(2))];['R^2 = ' num2str(R)]})
    
end