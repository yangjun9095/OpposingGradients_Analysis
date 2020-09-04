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
% projects = LoadMS2Sets('r3-new')

% alternative is using the prefixes = getProjectPrefixes(dataType,varargin) function
projectsForFit = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3',...
                    'r0_RuntNull','r1_RuntNull','r2_RuntNull','r3_RuntNull',...
                    'r1_mid_RuntNull','r1_close_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};


% Because Jonathan's script doesn't support tuning the time window of
% fitting, we need to check this pretty much one by one, to find a good
% time window of inference (For example, it really depends on the quality
% of the mean spot fluorescence traces, as well as how long we have taken
% the measurement, etc.)

typeindice = 1;
dataType = projectsForFit{typeindice}
projects = getProjectPrefixes(dataType,'onlyApproved')

%% for each dataset(Prefix), run the fitting script
for set = 1:length(projects)
    Prefix = cell2mat(projects(set));
    FitMeanVectorWholeCycleNC14('Prefix', Prefix, 'NCWindow', [3, 25])
end



%% For each DataType, check the fitting
% typeindice = 1;
dataType = projectsForFit{typeindice}
projects = getProjectPrefixes(dataType,'onlyApproved')

% for each dataset(Prefix), check the fitting, approve/disapprove the fit at each AP bin
for set = 1:length(projects)
    Prefix = cell2mat(projects(set));
    ApproveMeanVectorWholeCycleFitsNC14('Prefix', Prefix)
end

%% Extract the fitted values
