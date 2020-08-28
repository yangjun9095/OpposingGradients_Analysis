% projectName = 'OpposingGradients'
% Last updated : 8/30/2020
% Description : fitting a simple model of transcriptional cycle to our
% averaged MS2 traces to extract parameters as the initial slope, duration,
% etc.

% Inputs : 1) single embryo (Prefix), 2) DataType (LoadMS2Sets) or
% AveragedDatasets, 3) individual embryos (optional...)


%% Use Jonathan's script for fitting the whole NC14 with initial rise, steady-state, and exponential decay
% I need more description on which parameters are extracted from the model
% fitting.

% Grab Prefixes using LoadMS2sets 
% alternative is using the getPrefixesFromDataStatusTab.m function
data_r3 = LoadMS2Sets('r3-new')

%% for each dataset(Prefix), run the fitting script
for set = 1:length(data_r3)
    Prefix = data_r3(set).Prefix;
    FitMeanVectorWholeCycleNC14('Prefix', Prefix)
end

%% for each dataset(Prefix), check the fitting, approve/disapprove the fit at each AP bin
for set = 1:length(data_r3)
    Prefix = data_r3(set).Prefix;
    ApproveMeanVectorWholeCycleFitsNC14('Prefix', Prefix)
end

%% Extract the fitted values
