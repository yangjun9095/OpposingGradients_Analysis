% projectName = 'OpposingGradients'
% Last updated : 9/1/2020
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
                    'r1_close_RuntNull','r1_mid_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};


%% For each dataType, perfor the fitting process in the decay regime of NC14.                
% For each dataset, they might need their own time windows for fitting, as
% in [000] case.
for typeIndex = 1:length(DataTypesForFit)
    dataType = projectsForFit{typeIndex};
    projects = getProjectPrefixes(dataType,'onlyApproved');

    % for each dataset(Prefix), run the fitting script
    for set = 1:length(projects)
        Prefix = cell2mat(projects(set));
        fit_ExpDecay_NC14DecayRegime(Prefix,'TimeWindow',[0 20])
    end
end

% show the progress
%waitbar(typeindex/length(projectsForFit),['fitting in progress : ',num2str(typeindex),'/16 processed'])


%% For each DataType, try the fitting with numerically calculating the max*(1-1/e) point, rather than fitting.
% This is because some traces are dominated by the noise.

for typeIndex = 1:length(projectsForFit)
    dataType = projectsForFit{typeIndex};
    projects = getProjectPrefixes(dataType,'onlyApproved');

    % for each dataset(Prefix), run the fitting script
    for set = 1:length(projects)
        Prefix = cell2mat(projects(set));
        fit_ExpDecay_NC14DecayRegime_upto30min(Prefix)
    end
end

%% Extract the fitted values
