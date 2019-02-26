function ReplaceFieldsForCompiledNuclei
% Replacing the fields in CompiledNuclei (MeanVectorAP, and SDVectorAP) with values after the
% cytoplasmic fluo / Kg subtraction.
% This calculation is done in the previous step using the Runt_Cyto_Fluo.m
% script.

clear all
Prefix = '2018-05-30-Runt-JB3-MCP-mCherry-vasa-eGFP1';

% Load the original CompiledNuclei
load(['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\',Prefix,'\CompiledNuclei.mat'])
% Load the edited dataset (subtracted with background free FP inside
% nuclei)
Dataset = load(['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Runt_without_Cytoplasmic_Fluo\',...
            Prefix,'.mat'])

% Define MeanVectorAP and SDVectorAP
MeanVectorAP = Dataset.TF_Nuclei;
SDVectorAP = Dataset.TFError;

% Save all the variables for CompiledNuclei
save(['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\',...
        Prefix,'\CompiledNuclei.mat'],...
        'CompiledNuclei','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','ncFilterID','ncFilter','APbinID','APFilter',...
        'MeanVectorAP','SDVectorAP','NParticlesAP','MeanVectorAll',...
        'SDVectorAll','NParticlesAll','MaxFrame','MinAPIndex','MaxAPIndex',...
        'AllTracesVector','AllTracesAP',...
        'MeanCyto','SDCyto','MedianCyto','MaxCyto',...
        'MeanCytoAPProfile','SDCytoAPProfile','SECytoAPProfile', '-v7.3')
end