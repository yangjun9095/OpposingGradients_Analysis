function compile_InitialSlopes_AllConstructs
%% Description 
% Compiling and avearging initial slopes from multiple embryos
% Generate plots, and also save the fields that are calculated.

%This code is using the "mean fits" files generated by my version of FitMeanAPAsymetric 
%to plot the initial rates along APaxis. 
%The T_on (time of activation of the transcription,extrapolated from the initial rate)
%can also be plotted.

%% Make a structure that compiles all information
DataTypesForFit = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3',...
                    'r0_RuntNull','r1_RuntNull','r2_RuntNull','r3_RuntNull',...
                    'r1_close_RuntNull','r1_mid_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};

% Name of constructs so that we can label plots and files
constructNames = {'000','100','011','111','001','010','110','101',...
                    '000, null','100, null','011, null','111, null','001, null','010, null','110, null','101, null'};
                
% SAve the DataType and LoadMS2Sets into a structure, named compiledData
% The first row will contain the field names
compiledData{1,1} = 'DataType';
compiledData{1,2} = 'compiledData';

% Loop through all DataTypes to load comipled data from each dataType
for DataType = [1,9]%1:length(DataTypesForFit)
    compiledData{DataType+1,1} = DataTypesForFit{DataType};
    compiledData{DataType+1,2} = LoadMS2Sets(DataTypesForFit{DataType},'dontCompare');
end


%% Extract the fitted values from all DataTypes
% First, initialize the structure with the field names
compiledData{1,3} = 'fitted Rate';
compiledData{1,4} = 'fitted Rate_SD';
compiledData{1,5} = 'fitted T_ON';
compiledData{1,6} = 'T_peak';
compiledData{1,7} = 'Tau';
% compiledData{1,8} = 'Tau_SD';
compiledData{1,8} = 'Tau_numer';

% Go over the structure (compiledData) to extract the fitted initial slope,
% SD, number of embryos, T_on, T_peak, Tau, etc.

for i=1:length(DataTypesForFit)
    [fittedRate,fittedRateSD,fittedTon,T_peak, Tau, Tau_SD, Tau_numer] = Extract_Fields_MeanFits(compiledData{i+1,2},'Asymmetric');
    
    compiledData{i+1,3} = fittedRate;
    compiledData{i+1,4} = fittedRateSD;
    compiledData{i+1,5} = fittedTon;
    compiledData{i+1,6} = T_peak;
    compiledData{i+1,7} = Tau;
%     compiledData{i+1,8} = Tau_SD;
    compiledData{i+1,8} = Tau_numer; % numerically calculated from the maximum at 40min into nc14
end


%% Calculate the average using nanmean, nanstd

% Loop through all datasets to calculate the mean and SEM of the fitted
% rate, T_ON, Tau, and duration (in NC14)

% NOTE : the [000] and [000], Runt null are taken for less than 30 min,
% thus we need to use the 

% initialize the structure column
% fitted initial slope and SEM
compiledData{1,9} = 'fittedRate_mean';
compiledData{1,10} = 'fittedRate_SEM';
% fitted T_ON and SEM
compiledData{1,11} = 'fitted_T_ON_mean';
compiledData{1,12} = 'fitted_T_ON_SEM';
% Duration of Txn (T_peak - T_ON + Tau)
% average and SEM over multiple embryos of the same genotype
compiledData{1,13} = 'Duration_mean';
compiledData{1,14} = 'Duration_SEM';
% duration of individual embryos
compiledData{1,15} = 'Duration_individual';

NC = 3; %nc14

for i=1:length(DataTypesForFit)
    % clear the variables for each run
    vars = {'fittedRate_mean','fittedRate_SEM',...
                'fittedTon','fitted_T_ON_SEM','T_peak','Duration_individual',...
                'Duration_mean','Duration_SEM'};
    clear (vars{:})
    
    % initial slope
    fittedRate_mean = nanmean(compiledData{i+1,3}(:,NC,:),3); % only nc14
    fittedRate_SEM = nanstd(compiledData{i+1,3}(:,NC,:),0,3)./sqrt(length(compiledData{i+1,2})); 
    compiledData{i+1,9} = fittedRate_mean;
    compiledData{i+1,10} = fittedRate_SEM;
    
    % T_ON
    fitted_T_ON_mean = nanmean(compiledData{i+1,5}(:,NC,:),3); % only nc14
    fitted_T_ON_SEM = nanstd(compiledData{i+1,5}(:,NC,:),0,3)./sqrt(length(compiledData{i+1,2})); 
    compiledData{i+1,11} = fitted_T_ON_mean;
    compiledData{i+1,12} = fitted_T_ON_SEM;
    
    % Txn Duration 
    T_peak = compiledData{i+1,6};
    T_ON = squeeze(compiledData{i+1,5}(:,NC,:));
    Tau = compiledData{i+1,8}; % numerically calculated from the max @ 30min
    %Tau = compiledData{i+1,7}; % fitted : [000] and [000], Runt null (1st
    %and 9th datatype)
    
    % filter out the Tau value that is too long
    Tau(Tau>15) = nan; % threshold of 60 min, this should be revisited later.
    Tau(Tau<2) = nan; % threshold of 60 min, this should be revisited later.
    
    % Duration = (T_peak - T_ON + Tau)
    Duration_individual = T_peak - T_ON + Tau; % APbins x embryos
    Duration_mean = nanmean(Duration_individual,2);
    Duration_SEM = nanstd(Duration_individual,0,2)./sqrt(length(compiledData{i+1,2}));
    compiledData{i+1,13} = Duration_mean;
    compiledData{i+1,14} = Duration_SEM;
    compiledData{i+1,15} = Duration_individual;
end

%% check the duration of Txn over AP axis
% hold on
% % for i=1:5
% %     plot(0:0.025:1, Tau(:,i))
% % end
% errorbar(0:0.025:1, Duration_mean, Duration_SEM)
% xlabel('embryo length')
% ylabel('duration (min)')
% ylim([0 60])
% 
% % save the plot

%% Save the structure, compiledData for future usage (in plotting scripts)
save('S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\compiledData.mat',...
        'compiledData')

end