function RunCompileParticles
% This function is for running CompileParticles for multiple datasets
% easily. 

%% Input
% DataType in the tab of excel file, DataStatus.xls
% It uses LoadMS2Sets to get the Prefix info.
% Thus, we have to be careful about the ComputerFolders.csv
% that they don't have the sam
%% Output
% It runs CompileParticles, and saves CompiledParticles, etc. in the
% designated folder.

%% CompileParticles
% Define the dataset
DataType = 'r3';
Data = LoadMS2Sets(DataType);

for i=1:length(Data)
    Prefix = Data(i).Prefix;
    CompileParticles(Prefix,'SkipAll','ApproveAll','MinParticles',1)
end

%% AverageDatasets_FineTuning
% This is also based on LoadMS2Sets, thus needs to have only one
% DropboxFolder that has DataStatus.xls having the tab matches with the 'DataType'
DataType = 'r3';
Path = 'E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient';

AverageDatasets_FineTuning(DataType,'savePath',Path);

end
