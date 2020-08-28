%% Plotting Module

% The point is to put together blocks of codes to generate plots, as I've
% been writing the same code blocks over and over.

%% Color module
% color sets are important especially when there are >10 categories.
% 
%% Type1. 2D plot (errorbar/plot)

% define the figure handle, fig_name
fig_name = figure; 
hold on
errorbar(X, Y, Z, 'LineWidth',2,'Color',ColorChoice(1,:))

% xTicks, yTicks
xticks([0.2 0.3 0.4 0.5 0.6 0.7 0.8])

set(gca,'yticklabel',[])

% no title, no-caps on the axis labels
xlabel('')
ylabel('')

legend('','','Location','SouthWest')

box on

StandardFigure(fig_name, fig_name.CurrentAxes)

% Save the plot
figPath = '';
saveas(gcf,[figPath,filesep,'name','.tif']); 
saveas(gcf,[figPath,filesep,'name','.pdf']); 

%% Type 2. 3D plot (surface,scatter, etc.)

