function compare_txnOutput_WT_RuntNulls(DataType)
% Description
% This is a script to generate plots of comparing the Txn output
% with/without the repressor (Runt).

% Last updated : 5/1/2020
% Agenda : For a fair assessment of the Repressor's activity ( what it
% changes exactly: RNAP loading rates, duration, fraction ON, etc.), we
% need to compare one construct with/without repressor protein.

% For this, I'll divide the sections into several parts, depending on which
% metric that we're looking at.

% 1) instantaneous mean spot fluo
% 2) instantaneous RNAP loading rate (either inferred with Jonathan's MCMC
% script or his another three-phase fitting script).
% 3) instantaneous fraciton on
% 4) Txn duration time (Use Jonathan's script and refer to the Fig.S4 in
% Garcia 2013).

% Note that I need to repeat this for the [000]. But, from the raw traces,
% the fold-change looks like 1, for the instantaneous mean spot fluo,
% initial slope, Txn duration time. Fraction on should be checked further.


%% Load the datasets of [111] 
% Define the file path
DropboxPath = 'S:/YangJoon/Dropbox';
filePath = [DropboxPath,filesep,'OpposingGradient/OpposingGradients_ProcessedData/AveragedDatasets_Feb2020'];
Data_r3 = load([filePath, filesep, ])

end