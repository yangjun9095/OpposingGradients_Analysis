function supple_compare_differentEmbryos_sameConstruct(DataType)
% Description
% This script is for assessing the embryo-to-embryo variability.
% 
%% First, let's check whether all [000] new datasets male/female are equal.
% DataType names : 
% All r0-new construct (fixed from 
% r0-new : non-sexed embryos (5) up to nc14.
% r0-new-male : embryos only up to nc13 
% r0-new-female : embryos only up to nc13 
% r0-new-mixed : mix of males + female 

%% Load the processed datasets (by AverageDatasets)
% file path
filePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\TxnOutput_sexed';

% Load the datasets
Data_r0_new = load([filePath, filesep, 'r0-new.mat']);
Data_r0_new_male = load([filePath, filesep, 'r0-new.mat']);
Data_r0_new_female = load([filePath, filesep, 'r0-new.mat']);
Data_r0_new_mixed = load([filePath, filesep, 'r0-new.mat']);

%% First, let's check whether all [111] new datasets male/female, old vs new are equal.
% DataType names : 
% All r3 construct (old vs new)
% r3 : old (non-sexed)
% r3-new : 



%% Load the processed datasets (by AverageDatasets)
% file path
filePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\TxnOutput_sexed';

% Load the datasets
% Data_r0_new = load([filePath, filesep, 'r0-new.mat']);
% Data_r0_new_male = load([filePath, filesep, 'r0-new.mat']);
% Data_r0_new_female = load([filePath, filesep, 'r0-new.mat']);
% Data_r0_new_mixed = load([filePath, filesep, 'r0-new.mat']);
end