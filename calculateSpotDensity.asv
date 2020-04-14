function calculateSpotDensity(DataType,varargin)
% Description
% This function is to calculate the spot density 
% spot density = # of spots / area of APbins
% Caveats : 
% 1) we need to ignore APbins whose area is just too small, compared
% to the half of the maximum APbin area.
% 2) 

%% Method
% Get the number of spots per AP bin at each frame (nSpots)
% Get the APbin area for corresponding AP bins (areaAPbins)
% Then, get the ratio of "nSpots/areaAPbins"
% We call it as density. For specific datasets, we take the maximum density
% as 1 (this could be either maximum density across different datasets,
% etc.), assuming the fraction ON in the anterior (20%) is almost always 1.

% Another method is just going back to the raw movies, then do some
% curations of nuclear masks for a couple of frames in NC14 (or NC13 as
% well), then apply those masks for all of the frames.

%% Method1.

% Load the datasets
Data = LoadMS2Sets(DataType)


end