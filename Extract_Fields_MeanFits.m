%% Extract useful fields
function [fittedRate,fittedRateSD,fittedTon] = Extract_Fields_MeanFits(Data)
    % Data is the compiled datasets by LoadMS2Sets.m
    
    % From this, let's construct 3D matrices for fitted Rate and fitted T on,
    % which have dimensions like (AP,NC,index of embryo)

    % Initialize the matrices, fill with Nans.
    fittedRate = nan(41,3,length(Data));
    fittedRateSD = nan(41,3,length(Data));
    fittedTon = nan(41,3,length(Data));

    for i=1:length(Data)
        % check if there's MeanFitsV2 (which is from Asymmetric fit)
        if isfield(Data(i),'MeanFits')
            MeanFits = Data(i).MeanFits;
        elseif isfield(Data(i),'MeanFitsV2')
            MeanFits = Data(i).MeanFitsV2;
        end
        
        % Extract the fitted values.
        for j=1:length(MeanFits)
            for NC = 1:3 % NC12 to NC14 by default
                if ~isempty(MeanFits(j,NC).RateFit)&& (MeanFits(j,NC).Approved==1) % Only take the values from the ones exist, and also approved
                    fittedRate(j,NC,i) = MeanFits(j,NC).RateFit;
                    fittedRateSD(j,NC,i) = MeanFits(j,NC).SDRateFit;
                    fittedTon(j,NC,i) = MeanFits(j,NC).TimeStart;
                end
            end
        end
    end
end