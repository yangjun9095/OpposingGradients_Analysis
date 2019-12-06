function Supple_05_02_LlamaTagged_Runt_BGsubtracted
%% Test if the background subtraction has worked well.

%% Start with one dataset
Prefix = '2018-12-03-RuntJB3-vasa-eGFP-Pos8'

%% File Directory
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath,...
    configValues, movieDatabasePath, movieDatabaseFolder, movieDatabase]=...
    DetermineLocalFolders(Prefix);

%% Load the dataset
NucFluoData = load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'])

%% Extract useful fields for plotting
APaxis = 0:0.025:1;

for i=1:length(NucFluoData.ElapsedTime)
    clf
    hold on
    errorbar(APaxis, NucFluoData.MeanVectorAP(i,:), NucFluoData.SDVectorAP(i,:))
    errorbar(APaxis, NucFluoData.MeanVectorAP_BGsubtracted(i,:), NucFluoData.SDVectorAP_BGsubtracted(i,:))
    hold off
    pause 
end