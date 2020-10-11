function [Max, Min, DynamicRange, BoundaryPosition, Slope] = getGradientFeatures (X,Gradient)
% DESCRIPTION
% This function grabs inputs of a gradient over an AP axis (or could be
% DV), then calculates 
% 1) Maximum, Minimum, and dynamic range (Max - Min), 
% 2) Boundary position (half-maximum)
% 3) Slope (which is a slope of a tangential line from the boundary
% position)

% OUTPUT
% 1) Maximum, Minimum, and dynamic range (Max - Min), 
% 2) Boundary position (half-maximum)

% INPUT
% X : X-axis (either AP axis or DV)
% Gradient : Gradient of expression, could be either protein or mRNA, or
% P_bound, etc.

%% Define different outputs (Max, Min, and Dynamic range)
Max = max(Gradient);
Min = min(Gradient);
DynamicRange = Max - Min;

%% Get boundary position, and slope
MidValue = Min + 0.5*DynamicRange;
Residuals = (Gradient - MidValue).^2;
SlopeTemp = [0 diff(Gradient)];
slopeSign = (SlopeTemp < 0);
%slopeSign(slopeSign ==0) = nan; % Make elements of slopeSign (positive slopes) as nans, for easier calculation.

% Residual x Slope = estimate of distance from the half-maximum (y-axis),
% so finding the minimum would give us the boundary position.
Residual_with_sign = Residuals.*slopeSign;
Residual_with_sign(Residual_with_sign==0) = nan;

MinResidual = min(Residual_with_sign);
BoundaryIndex = find(Residual_with_sign==MinResidual); 
if ~isnan(BoundaryIndex)
    BoundaryIndex = BoundaryIndex(1); % Pick the first of that boundary
    BoundaryPosition = X(BoundaryIndex);
else 
    BoundaryPosition = nan;
end

if isempty(BoundaryPosition)
    BoundaryPosition = nan;
end

%Medium = min((Gradient - MidValue).^2);
%BoundaryIndex = find(Residuals == MidTemp);
%BoundaryPosition = X(BoundaryIndex);

% Get a tangential line from the mid-point
spacing = 2; % this is to define how much of x-window to do a linear fit.
window = BoundaryIndex-spacing:BoundaryIndex+2;
% In case the fitting doesn't work
try
    fittedResult= polyfit(X(window),Gradient(window),1);
    Slope = fittedResult(1);
catch
    Slope = nan;
end

end