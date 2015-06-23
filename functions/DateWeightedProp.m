function mean_optprop = DateWeightedProp(analysis_date, optprops)
    unique_date = unique(analysis_date(~isnan(optprops)));
    analysis_weight = zeros(size(analysis_date));
    for R1 = 1:length(unique_date)
        n = length(find((analysis_date==unique_date(R1).*~isnan(optprops))));
        analysis_weight(analysis_date==unique_date(R1)) = 1/n;
    end
    mean_optprop = nansum(analysis_weight.*optprops)/length(unique_date);
end