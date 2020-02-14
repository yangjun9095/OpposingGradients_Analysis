function modeling_01_generatePredictions_hbP2
% DESCRIPTION
% Script for generating predictions for 2 Runt sites or more, 
% based on 0, 1 binding site configurations

% The goal here is to trying to generate predictions of hb P2 + 2 Runt
% sites or more based on 0, and 1 sites (at different positions).

% For this, we don't know how this system works actually, meaning the way
% repressor works, also it doesn't seem like that simple (not a direct
% competition, nor a direct repression to the RNAP recruitment). But, from
% our theoretical exploration, it doesn't seem to be distinguishable by
% modeling. 
% Thus, maybe it's better to make this script to be able to compare the
% results from different models 
% : not just for adding complexity (from Hill model to individual sites), 
% but also deal with different scenarios of repression (not sure if this is necessary)

%% Step0. Can we do somewhat phenomenological Hill model like in Jeehae's paper?
% I don't know how to incorporate the repressor terms in here...
% How about 1/(1 + ([R]/Kr)^Nr)

%% Step1. Let's start with the simplest Hill model of Bcd activating the hb P2

end