function Useful_Modules_collection
% Collection of useful modules

%% Module 1. Color gradation plot
% Indexing -  jStart:jEnd

% Different types of color maps
% From Meghan, we now have more linear, and color-blind friendly color
% maps, such as viridis, inferno, magma, and plasma, etc. They are all in
% our mRNADynamics folder, and we can call them easily.

jStart = 1;
jEnd = 255;
colormap(jet(256)); % or viridis, inferno, magma, or plasma
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for i=iStart:iEnd
    plot(X,Y,'color',Color(i-jStart+1,:))  %need to consider background level
end

%% Module 2. Color assignments as PBoC
%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;AzzzzZXXxxxxxx
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179.155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

ColorChoice = [colorDict.blue; colorDict.red; colorDict.green; colorDict.purple]; % 4 embryos max. it could be extended easily
lineColor = ['b', 'r', 'g', 'p'];
end