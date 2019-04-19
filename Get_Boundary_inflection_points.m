function [BoundaryPosition, Slope] = Get_Boundary_inflection_points(X,Gradient)
% DESCRIPTION
% This function is to get boundary features for a "Gradient" over "X" as
% Tao did in Chen, 2012.
% The definition is to find inflection points over all X axis, then get
% the median for the "Boundary position".

%% Get the double derivative
window_smooth = 5;
Gradient_smoothend = movmean(Gradient,window_smooth);

single_derivative = diff(Gradient_smoothend);
double_derivative_Gradient = diff(diff(Gradient_smoothend));

% Question 1 : Do we need interpolation?

% Question 2 : Should we use smoothened data or raw? I'd use Smoothened

% Questuib 3 : Is the boundary position (like the positions of inflection
% points affected by the averaging window ?
% Define the window of inflection points that we need to find, 
% I'll define it as 20~60% for now, since the inflection points in anterior
% doesn't mean anything.
%Window = 1:length(X);
inflection_points = find(double_derivative_Gradient<0.0001);

%% Plot to check
% hold on
% plot(X, Gradient)
% plot(X, Gradient_smoothend)
% plot(X(inflection_points),Gradient_smoothend(inflection_points),'or','MarkerSize', 10)

%% find Median

BoundaryPosition = X(floor(median(inflection_points)));
Slope = single_derivative(floor(median(inflection_points)))/0.01;
end
