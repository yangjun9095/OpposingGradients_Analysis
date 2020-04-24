function calculateSpotDensity(DataType,varargin)
% Description
% This function is to calculate the spot density 
% spot density = # of spots / area of APbins
% Caveats : 
% 1) We need to ignore APbins whose area is just too small, compared
% to the half of the maximum APbin area.
% 2) 
% Assumptions (that could be relaxed later)
% 1) We will only think about NC14 as our current Runt null datasets have
% only have NC14 intact. Thus, we will synchronize all the vectors to the
% beginning of NC14 (In other words, the CheckDivisionTimes.m is very important).
% 2) For a long enough time series, we will assume 

%% Method
% Get the number of spots/particles per AP bin at each frame (nParticles)
% Get the APbin area for corresponding AP bins (areaAPbins)
% Then, get the ratio of "nSpots/areaAPbins"
% We call it as density. For specific datasets, we take the maximum density
% as 1 (this could be either maximum density across different datasets,
% etc.), assuming the fraction ON in the anterior (20%) is almost always 1.

% Another method is just going back to the raw movies, then do some
% curations of nuclear masks for a couple of frames in NC14 (or NC13 as
% well), then apply those masks for all of the frames.

%% Method1.

% Temporary
DataType = 'r3_RuntNull';

% Load the datasets
Data = LoadMS2Sets(DataType)

% Extract useful fields
NParticles = Data.NParticlesAP;

end