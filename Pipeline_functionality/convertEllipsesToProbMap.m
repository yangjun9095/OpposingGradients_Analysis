% Convert Ellipses To Probability map
function convertEllipsesToProbMap(Prefix)
% This script is for the case when we have a nice segmentation of nuclei,
% but bad tracking. So, we will use the Ellipses.mat to generate a fake
% probability map, which can be fed into the Tr2D for a better tracking.

%% Load the dataset
% We need nuclei segmentation(Ellipses.mat), pixel size, radius, etc.
% 

%% Make a binary image of nuclei/cytoplasm

%% Save

end
