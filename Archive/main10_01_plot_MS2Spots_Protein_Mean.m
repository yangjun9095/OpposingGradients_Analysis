function main10_01_plot_MS2Spots_Protein_Mean(DataType, varagin)
%% DESCRIPTION

%% Options
% Use varargin for making options to be fed into this.
% i.e. NC, number of embryos, etc.

%% Load datasets
DataType = 'r3-RuntJB3-MCP-mCherry'; % as an example, this should be an input for this function!

Data = LoadMS2Sets(DataType);

%% Loop over multiple embryos

%% For a single embryo
% 1) extract useful fields
Particles = Data.Particles;
Nuclei = Data.Nuclei;

MesnSpotFluo = Particles.MeanVectorAP{1,1};
SDSpotFluo = Particles.SDVectorAP{1,1};
NSpots = Particles.NParticlesAP{1,1};


MeanNucFluo = Nuclei.MeanVectorAP;
SDNucFluo = Nuclei.SDVectorAP;
NNuclei = Nuclei.NParticlesAP;


end