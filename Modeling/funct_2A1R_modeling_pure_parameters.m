function funct_2A1R_modeling_pure_parameters
% This is for modeling 2A1R case, with pure sets of parameters


% Parameters :
% a : [A]/Ka
% r : [R]/Kr
% p = [P]/Kp
% w_ap : Activator on RNAP
% w_aa : Activator cooperativity
% w_rp : Rpressor on RNAP

%% Define the parameters
% The scheme is to set the parameters in certain regime, then fix all
% parameters except one, then tune one parameter to see how much it
% changes...
a = 10.^[(-3) (-2) (-1) 0 1 2 3];
r = 10.^linspace(-3,3,100);%[(-3) (-2) (-1) 0 1 2 3];
p = 10.^[(-3) (-2) (-1) 0 1 2 3];

w_ap = 10.^[0 1 2 3 ];
w_rp = 10.^[0 -1 -2 -3];
w_aa = 10.^[0 1 2 3];

% a = a(6);
% r = r;
% p = p(3);
% 
% w_ap = w_ap(2);
% w_rp = w_rp(2);
% w_aa = w_aa(2);

for i=1:length(a)
    for j=1:length(p)
        for k = 1:length(w_ap)
            for l = 1:length(w_rp)
                for m = 1:length(w_aa)
                    clear Numerator
                    clear Denominator
                    Numerator = a(i).^2.*p(j).*w_ap(k).^2.*w_aa(m).*(1 + r.*w_rp(l));
                    Denominator = (1+r).*(1+a(i).^2.*w_aa(m)) + a(i).^2.*p(j).*w_ap(k).^2.*w_aa(m).*(1 + r.*w_rp(l));

                    Pbound{i,j,k,l,m} = Numerator./Denominator;
                end
            end
        end
    end
end

%% Change one parameter at one time, fix everything else.
% First, let's tune activator conc. (a)
% Fix p, w_ap, w_rp, w_aa
hold on
for i=1:length(a)
    semilogx(r,Pbound{i,2,2,2,2})
end

set(gca,'Xscale','log')
set(gca,'Yscale','log')

title('P_{bound} vs Repressor (for various [A]/K_{A})')
xlabel('[Repressor]/K_{R}')
ylabel('P_{bound}')
legend('10^{-3}','10^{-2}','10^{-1}','1','10','10^{2}','10^{3}')

StandardFigure(gcf,gca)

% Save the figure
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_2A1R_ParameterExploration';
saveas(gcf,[FigPath,filesep,'P_bound_2A1R','p=',num2str(p(3)),'_w_ap=',num2str(w_ap(2)),'_w_rp=',num2str(w_rp(2)),'_w_aa=',num2str(w_aa(2)) , '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_2A1R','p=',num2str(p(3)),'_w_ap=',num2str(w_ap(2)),'_w_rp=',num2str(w_rp(2)),'_w_aa=',num2str(w_aa(2)) , '.pdf']); 
%% Second, tuning the w_ap
% Fix a, w_ap, w_rp, w_aa
% a = 10;
% p = 0.01;
% w_rp = 0.1;
% w_aa = 10;

hold on
for k=1:length(w_ap)
    semilogx(r,Pbound{5,2,k,2,2})
end

set(gca,'Xscale','log')
set(gca,'Yscale','log')

title('P_{bound} vs Repressor (for various w_{ap})')
xlabel('[Repressor]/K_{R}')
ylabel('P_{bound}')
legend('1','10','10^{2}','10^{3}')
ylim([0 1.2])
StandardFigure(gcf,gca)

% Save the figure
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_2A1R_ParameterExploration';
saveas(gcf,[FigPath,filesep,'P_bound_2A1R_tuning_w_ap','_a=',num2str(a(5)),'p=',num2str(p(3)),'_w_rp=',num2str(w_rp(2)),'_w_aa=',num2str(w_aa(2)) , '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_2A1R_tuning_w_ap','_a=',num2str(a(5)),'p=',num2str(p(3)),'_w_rp=',num2str(w_rp(2)),'_w_aa=',num2str(w_aa(2)) , '.pdf']); 

%% Third, tuning the w_rp
% Fix a, w_ap, w_rp, w_aa
% a = 10;
% p = 0.01;
% w_ap = 10;
% w_aa = 10;

hold on
for l=1:length(w_rp)
    semilogx(r,Pbound{5,2,2,l,2})
end

set(gca,'Xscale','log')
set(gca,'Yscale','log')

title('P_{bound} vs Repressor (for various w_{rp})')
xlabel('[Repressor]/K_{R}')
ylabel('P_{bound}')
legend('1','10^{-1}','10^{-2}','10^{-3}')

ylim([0 1.2])
StandardFigure(gcf,gca)

% Save the figure
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_2A1R_ParameterExploration';
saveas(gcf,[FigPath,filesep,'P_bound_2A1R_tuning_w_rp','_a=',num2str(a(5)),'p=',num2str(p(3)),'_w_ap=',num2str(w_ap(2)),'_w_aa=',num2str(w_aa(2)) , '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_2A1R_tuning_w_rp','_a=',num2str(a(5)),'p=',num2str(p(3)),'_w_ap=',num2str(w_ap(2)),'_w_aa=',num2str(w_aa(2)) , '.pdf']);

%% Fourth, tuning the w_aa
% Fix a, w_ap, w_rp, w_aa
% a = 10;
% p = 0.01;
% w_ap = 10;
% w_rp = 0.1;

hold on
for m=1:length(w_aa)
    semilogx(r,Pbound{5,2,2,2,m})
end

set(gca,'Xscale','log')
set(gca,'Yscale','log')
title('P_{bound} vs Repressor (for various w_{aa})')
xlabel('[Repressor]/K_{R}')
ylabel('P_{bound}')
legend('1','10^{1}','10^{2}','10^{3}')

ylim([0 1.2])
StandardFigure(gcf,gca)

% Save the figure
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_2A1R_ParameterExploration';
saveas(gcf,[FigPath,filesep,'P_bound_2A1R_tuning_w_aa','_a=',num2str(a(5)),'p=',num2str(p(3)),'_w_ap=',num2str(w_ap(2)),'_w_rp=',num2str(w_rp(2)) , '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_2A1R_tuning_w_aa','_a=',num2str(a(5)),'p=',num2str(p(3)),'_w_ap=',num2str(w_ap(2)),'_w_rp=',num2str(w_rp(2)) , '.pdf']);

%% Fifth, tuning the p
% Fix a, w_ap, w_rp, w_aa
% a = 10;
% p = ?
% w_ap = 10;
% w_rp = 0.1;
% w_aa = 10;

hold on
for j=1:length(p)
    semilogx(r,Pbound{5,j,2,2,2})
end

set(gca,'Xscale','log')
set(gca,'Yscale','log')

title('P_{bound} vs Repressor (for various p)')
xlabel('[Repressor]/K_{R}')
ylabel('P_{bound}')
legend('10^{-3}','10^{-2}','10^{-1}','1','10^{1}','10^{2}','10^{3}')

ylim([0 1.2])
StandardFigure(gcf,gca)

% Save the figure
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_2A1R_ParameterExploration';
saveas(gcf,[FigPath,filesep,'P_bound_2A1R_tuning_p','_a=',num2str(a(5)),'p=',num2str(p(3)),'_w_ap=',num2str(w_ap(2)),'_w_rp=',num2str(w_rp(2)), '_w_aa=',num2str(w_aa(2)), '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_2A1R_tuning_p','_a=',num2str(a(5)),'p=',num2str(p(3)),'_w_ap=',num2str(w_ap(2)),'_w_rp=',num2str(w_rp(2)), '_w_aa=',num2str(w_aa(2)), '.pdf']); 

%% Part2. Fold-Change
% I'll use the derivation that I put on the Overleaf document for the
% fold-change (2 activator sites + 1 repressor site)
% From the expression calculated from Overleaf document, with many
% assumptions.
% 2A1R case
for i=1:length(w_rp)
    FC_2A1R(i,:) = (1 + r*w_rp(i))./(1+r);
end

%% Plot the fold-change (r)
hold on
for i=1:length(w_rp)
    semilogx(r,FC_2A1R(i,:))
end

set(gca,'Xscale','log')
set(gca,'Yscale','log')

title('Fold-change vs Repressor')
xlabel('[Repressor]/K_{R}')
ylabel('Fold_Change')
legend('1','10^{-1}','10^{-2}','10^{-3}')

ylim([0 1.2])
StandardFigure(gcf,gca)

% Save the figure
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_2A1R_ParameterExploration';
saveas(gcf,[FigPath,filesep,'Fold_change_Repressor_w_rp','.tif']); 
saveas(gcf,[FigPath,filesep,'Fold_change_Repressor_w_rp','.pdf']); 

end