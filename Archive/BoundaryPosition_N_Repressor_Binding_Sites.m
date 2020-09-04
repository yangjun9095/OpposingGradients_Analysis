function BoundaryPosition_N_Repressor_Binding_Sites
% Generate the boundary position vs N_{R} (number of repressor binding
% sites)

% N_R
N_R = [0 1 2 3 ];
BoundaryPosition = [0.39 0.38 0.34 0.29];
BoundaryPosition_Chen_2012 = 1- [0.549 0.576 0.657 0.695];

BoundaryPosition_data_figure = figure;
hold on
plot(N_R,BoundaryPosition,'o')
plot(N_R,BoundaryPosition_Chen_2012,'o')
title('Boundary Position for different N_{R}')
xlabel('N_{R} : Number of repressor binding sites')
ylabel('Boundary position (Embryo Length)')
legend('Our measurement','Chen 2012')
xlim([0 7])
ylim([0 0.6])
standardizeFigure(gca,legend,[])

%% save the plot
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Figures-OpposingGradients\';
saveas(BoundaryPosition_data_figure,[FigPath 'Boundary_position_for_N_R' , '.pdf']); 
end