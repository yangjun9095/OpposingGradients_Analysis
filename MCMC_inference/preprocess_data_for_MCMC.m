function preprocess_data_for_MCMC
%% Description
% This script is for processing the data structure (compiledData.mat) into
% another structure that can be fed into the MCMC inference protocol.
% There could be many ways to do this, 
% 1) save WT data and Runt null data separately, with AP bins in another
% field.
% So, it'd be the "data" as the main structure, with fields like
% constructName, APbins, Rate_null, Rate_WT, Rate_null_individual, Rate_WT_individual,
% etc.

%% Load the compiledData.mat
FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';
load([FilePath, filesep, 'compiledData.mat'])

% Initialize the structure for the pre-processed data
data.constructName = {};
data.APbins = {};
% rates (averaged over embryos)
data.Rate_WT = {};
data.Rate_null = {};

% rates (individual embryos)
data.Rate_WT_individual = {};
data.Rate_null_individual = {};

%% Process one dataset as an example : 
% [001] (r1-close) : 5th element
for index = 1:8
    data(index).constructName = compiledData{index+1,1};
    data(index).APbins = [0:0.025:1];

    data(index).Rate_WT = compiledData{index+1,9}; % 9th column for the initial rate averaged over embryos
    data(index).Rate_null = compiledData{index+1+8,9};

    % for the initial slope from the individual embryos, we will only take the
    % NC14, which is the 3rd column from numAPbins x NC x numEmbryos
    data(index).Rate_WT_individual = squeeze(compiledData{index+1,3}(:,3,:));
    data(index).Rate_null_individual = squeeze(compiledData{index+1+8,3}(:,3,:));
end

%% save the structure "data"
save([FilePath, filesep,'PreProcessedData_ForMCMC.mat'],...
        'data')
end