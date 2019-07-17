% Script for re-running segmentSpots with former way of calculating the
% offset (flat plane).
% This script can be used for re-running the CompileParticles for 3D
% gaussian fitted intensity.

% I need to re-run 1) segmentSpots, 2) TrackmRNADynamics,
% 3)AddParticlePosition, 4) CompileParticles
% fit3DGaussianToAllSpots (optional).

%% Load the datasets.
% For each DataType, I'll use the LoadMS2Sets to load the name of the
% datasets.
DataTypes = {'r0-new-female'}%'r3-new-female','r3-new-male'}%{'r0-new-female','r1-new-female','r2-new-female','r3-new-female'...
                %'r0-new-male','r1-new-male','r2-new-male','r3-new-male'};%

%%

for i=1:length(DataTypes)
    clear Data
    clear f
    
    DataType = DataTypes{i};
    Data = LoadMS2Sets(DataType,'dontCompare');
    
    %length(Data);
        
    % Waitbar
    %f = waitbar(0,[DataType,' being analyzed']);
    % Get the SetName using indexing. Define the Prefix here.
    
    for j=2:length(Data)
        clear Prefix
        Prefix = Data(j).SetName(11:end-1)
        
        %waitbar(1/length(Data)*j,f,[Prefix, ' being anlayzed']);

        % Running segmentSpots
        % These options are only for specific datasets from the Opposing
        % gradients, should be adjusted as an option in the future.
        segmentSpots(Prefix,1040,'Shadows',1,'keepProcessedData')
        
        % TrackmRNADynamics
        TrackmRNADynamics(Prefix,'noRetracking')
        
        % AddParticlePosition
        AddParticlePosition(Prefix,'yToManualAlignmentPrompt')
        
        % CompileParticles
        CompileParticles(Prefix,'SkipAll','ApproveAll','MinParticles',1,'yToManualAlignmentPrompt')
    end
end   
    
    
    
    
    
    
    
    
    
    