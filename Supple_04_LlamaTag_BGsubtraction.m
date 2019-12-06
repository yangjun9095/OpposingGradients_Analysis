function Supple_04_LlamaTag_BGsubtraction
% Script to process raw LlamaTagged TF fluorescence to extract the TF-eGFP
% population.
% Here, I'll try different approaches to subtract free eGFP population, and
% see how rigorous we should be.

%% Approach 1. Subtracting nuc fluo of very anterior/posterior bins
% from nuc fluo, assuming that there's no TF in those regions. This requires
% previous knowledge about spatial gene expression patterns, for example,
% where the gene is expressed or not.

%% Load the datasets (synchronized by AverageDatasets_NuclearProtein.m)
% File path for calling the datasets
FilePath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

Runt_female = load([FilePath, filesep, 'Runt-1min-200Hz-Female-Averaged.mat']);

Runt_male = load([FilePath, filesep, 'Runt-1min-200Hz-Male-Averaged.mat']);

NoNB = load([FilePath, filesep, 'NoNB-Averaged.mat']);

%% Extract useful fields

% Female
MeanFluo_female = Runt_female.MeanVectorAP;
SDFluo_female = Runt_female.SDVectorAP;
NNuclei_female = Runt_female.NParticlesAP;
SEFluo_female = SDFluo_female./sqrt(NNuclei_female); % for SEM of individual embryos.
Time_female = Runt_female.ElapsedTime;
NC13_female = Runt_female.nc13;
NC14_female = Runt_female.nc14;
end