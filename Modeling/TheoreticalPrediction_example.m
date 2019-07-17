% 

X = [0 1 2 3 ];

Y = [10 5 4 3];

plot (X, Y, 'o')

title('Prediction of transcription rate')
xlabel('Number of repressor binding sites')
ylabel('Rate of transcription')
xlim([-1 4])
ylim([0 12])
xticks([0 1 2 3 4])

FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\';
%standardizeFigure_YJK(gca,legend)
saveas(gcf,[FigPath 'Theoretical_Prediction_schematics'  , '.pdf']); 



