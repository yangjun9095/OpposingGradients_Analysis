function PlotBleachingTest_MCP_mCherry
% This is for finding the imaging condition that does not bleach the
% MCP-mCherry. The idea is to use sequential scanning and ROI to have one
% embryo imaging with two different frame rates.

% Load the dataset
MCPmCherry = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-09-30-P2PV1-2xMCP-mCherry-ROI-BleachingTest-800Hz-400Hz-BetweenLines\CompiledParticles.mat')

%% Plot ROI vs non-ROI

% AP bin
AP = 11;

hold on
errorbar(MCPmCherry.ElapsedTime,MCPmCherry.MeanVectorAP_ROI(:,AP),MCPmCherry.SDVectorAP_ROI(:,AP))
errorbar(MCPmCherry.ElapsedTime,MCPmCherry.MeanVectorAP_nonROI(:,AP),MCPmCherry.SDVectorAP_nonROI(:,AP))

title('bleaching test')
xlabel('Time (min)')
ylabel('Mean spot fluorescence')
legend('ROI(800Hz)','nonROI(400Hz)')



end