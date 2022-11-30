% Compare different methods to subtract the LlamaTag background

% DESCRIPTION
% there are multiple ways to treat the LlamaTag background fluo,
% first is using the No nanobody datasets as free eGFP/mCherry
% second is using the cytoplasmic fluo as described in the Bothma, 2018

%% 1. No Nanobody method
FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';
Runt_female = load([FilePath, filesep, 'Runt_processed_female.mat']);
Runt_male = load([FilePath, filesep, 'Runt_processed_male.mat']);

%% Time-average for the 0-10 min into nc14
tLength = length(Time_female) - NC14_female + 1;
nAPbins = 41;
numEmbryos = numEmbryos_female + numEmbryos_male;
MeanFluo_mixedSex = nan(tLength, nAPbins, numEmbryos);

MeanFluo_mixedSex(:,:,1:4) = MeanFluo_female_BGsubtracted(NC14_female:end,:,:);

MeanFluo_mixedSex(:,:,5:8) = MeanFluo_male_BGsubtracted(NC14_male:end,:,:);
%% Plot the spatial gradient

%% 2. Cyto fluo method
% This is already done using the scripts deposited in the main repository.
% The result for the mixed sex is saved in 