function Useful_Modules_collection
% Collection of useful modules

%% Module 1. Color gradation plot
% Indexing -  jStart:jEnd

% Different types of color maps
% From Meghan, we now have more linear, and color-blind friendly color
% maps, such as viridis, inferno, magma, and plasma, etc. They are all in
% our mRNADynamics folder, and we can call them easily.
colormap(jet(256)); % or viridis, inferno, magma, or plasma
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for i=iStart:iEnd
    plot(X,Y,'color',Color(i-jStart+1,:))  %need to consider background level
end
end