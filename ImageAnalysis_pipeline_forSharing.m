% Sequence of pipeline scripts to run for the image analysis of OpposingGradients datasets
% Last updated : 5/1/2020 by Yang Joon Kim
%

%% Export the LIF file
% Note that now we don't have to export the .lif files to .tif files with
% LASX.
Prefix = ExportDataForLivemRNA
% Click the folder containing the raw movie (.lif file, .lif file does not
% have to be exported.)

%% filtering movie
filterMovie(Prefix,'highPrecision','customFilter','Difference_of_Gaussian_3D','nWorkers',1)

%% spot segmentation
% Thresh = ? % after ImageJ
segmentSpots(Prefix, Thresh,'Shadows',1,'keepProcessedData')%'fit3D','displayFigures')

%% TrackNuclei

%% TrackmRNADynamics

%% CompileParticle
CompileParticles(Prefix,'skipAll','ApproveAll','MinParticles',2,'MinTime',2)