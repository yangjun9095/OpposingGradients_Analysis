% Modeling / Simulation for justifying Hill function expression of Bcd
% cooperativity
% Last edited : Sept.2019
% YJK
function BuildingUp_Repression_Aztec_Pyramid_StatMech_model

%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

ColorChoice = [colorDict.blue; colorDict.red; colorDict.green; colorDict.purple;...
                colorDict.brown; colorDict.darkgreen; colorDict.yellow; colorDict.cyan]; % 4 embryos max. it could be extended easily
%lineColor = ['b', 'r', 'g', 'p'];
%% 1. 1A1R, R=0
%% Parameter sensitivity analysis
% The scheme will be picking up random set of parameters for initial
% points, like 100 times or so,
% then try to find the best fit using lsqnonlin. Then, see whether the
% outcome converges.

% First, let's define the profile that we're going to fit.
% p = 0.01
% w_ap = 10
% w_rp = 0.01

% P_bound(p, a, r, w_ap, w_rp)

Profile=P_bound_1A1R_NoRepressor([0.01, 10, 1]);

% semilogx(a,Profile)

% Second, initialize the values for lsqnonlin by randomly picking set of
% parameters.
% Set the range of parameters to be drawn
p = linspace(10^(-3),10^(3),10^5);
w_ap = linspace(10^(-3),10^(3),10^5);
Ka = linspace(10^(-1), 10, 10^5);

% Initialize the parameters
N = 100; % # of repeats
s = RandStream('mlfg6331_64'); % for randomization?

p0 = datasample(s,p,100);
w_ap0 = datasample(s,w_ap,100);
Ka0 = datasample(s,Ka,100);

%% Repeat for N times for finding the best fit
lb = [min(p), min(w_ap), min(Ka) ];
ub = [max(p), max(w_ap), max(Ka) ];
options = optimset('Display','iter');

for i=1:N
    x0 = [p0(i), w_ap0(i), Ka0(i)];
    fun= @(x)P_bound_1A1R_NoRepressor(x) - Profile;

    x(i,:) = lsqnonlin(fun, x0, lb, ub, options);
end

%% Plot to check the convergence
hold on
scatter3(x(:,1),x(:,2),x(:,3),'ob')
scatter3(0.01,10,1,100,'or')

% set(gca,'XScale','log')
% set(gca,'YScale','log')
% set(gca,'ZScale','log')
xlim([0 0.02])
ylim([0 11])
zlim([0 1.1])

title('Parameter sensitivity test')
xlabel('P')
ylabel('w_{ap}')
zlabel('K_{a}')
grid on
view(10,20)

% Save the figures
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_1A1R_ParameterSensitivity';
saveas(gcf,[FigPath,filesep,'P_bound_','parameter_sensitivity' , '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_','parameter_sensitivity' , '.pdf']);

%% Example fitted result and original profile
hold on
plot(a,Profile,'color',ColorChoice(1,:),'LineWidth',10)
for i=1:100
plot(a,P_bound_1A1R_NoRepressor(x(i,:)),'color',ColorChoice(2,:))
end
set(gca,'XScale','log')

title('P_{bound} profile')
xlabel('Activator (AU)')
ylabel('P_{bound}')
legend('Original','Fit','Location','northwest')
StandardFigure(gcf,gca)

% Save figures
saveas(gcf,[FigPath,filesep,'P_bound_','parameter_sensitivity_exampleFits' , '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_','parameter_sensitivity_exampleFits' , '.pdf']); 

%% 2. 1A1R, Repressor present.
%% Parameter sensitivity analysis
% The scheme will be picking up random set of parameters for initial
% points, like 100 times or so,
% then try to find the best fit using lsqnonlin. Then, see whether the
% outcome converges.

% First, let's define the profile that we're going to fit.
% p = 0.01
% w_ap = 10
% r = 10 (or 1?)
% w_rp = 0.01

% P_bound_1A1R(p, w_ap, Ka, r, w_rp)
InitialValue_1A1R = [0.01, 10, 1, 1, 0.1];

Profile=P_bound_1A1R([0.01, 10, 1, 1, 0.1]);

semilogx(a,Profile)

% Second, initialize the values for lsqnonlin by randomly picking set of
% parameters.
% Set the range of parameters to be drawn
p = linspace(10^(-3),10^(3),10^5);
w_ap = linspace(10^(-3),10^(3),10^5);
Ka = linspace(10^(-1), 10, 10^5);
r = linspace(10^(-2),10^2, 10^4);
w_rp = linspace(10^(-3),10^(3),10^5);

% Initialize the parameters
N = 100; % # of repeats
s = RandStream('mlfg6331_64'); % for randomization?

p0 = datasample(s,p,100);
w_ap0 = datasample(s,w_ap,100);
Ka0 = datasample(s,Ka,100);
r0 = datasample(s,r,100);
w_rp0 = datasample(s,w_rp,100);

%% Repeat for N times for finding the best fit
lb = [min(p), min(w_ap), min(Ka), min(r), min(w_rp) ];
ub = [max(p), max(w_ap), max(Ka), max(r), max(w_rp) ];
options = optimset('Display','iter');

for i=1:N
    y0 = [p0(i), w_ap0(i), Ka0(i), r0(i), w_rp0(i)];
    fun= @(y)P_bound_1A1R(y) - Profile;

    y(i,:) = lsqnonlin(fun, y0, lb, ub, options);
end

%% Plot to check the convergence
% This time, the number of varialbes are 5, so it's hard to plot...
% Let's calculate the distance between two points (orignal vs fitted) in
% 5-D space, then plot that distance.
for i=1:N
    Distance(i) = sqrt(sum((InitialValue_1A1R - y(i,:)).^2));
end

hold on
%scatter3(y(:,1),y(:,2),Distance,'ob')
scatter3(y(:,1),y(:,2),y(:,3),'ob')
% scatter3(InitialValue_1A1R(4), InitialValue_1A1R(5), 0, 100,'or')
scatter3(InitialValue_1A1R(1), InitialValue_1A1R(2), InitialValue_1A1R(3), 100,'or')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'ZScale','log')
% xlim([0 0.02])
% ylim([0 11])
% zlim([0 1.1])
% 
title('Parameter sensitivity test')
xlabel('P')
ylabel('w_{ap}')
zlabel('K_{a}')
grid on
view(10,20)
StandardFigure(gcf,gca)
% Save the figures
% FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_1A1R_ParameterSensitivity';
saveas(gcf,[FigPath,filesep,'P_bound_1A1R_RepressorPresent_','parameter_sensitivity' , '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_1A1R_RepressorPresent_','parameter_sensitivity' , '.pdf']); 

%% Continued.
% What is I fix the parameters regarding activators, and only tune the
% repressor parameters? Would it constrain the repressor parameters, r,
% w_{rp} better?


%% Example fitted result and original profile
hold on
plot(a,Profile,'color',ColorChoice(1,:),'LineWidth',10)
for i=1:100
plot(a,P_bound_1A1R_NoRepressor(x(i,:)),'color',ColorChoice(2,:))
end
set(gca,'XScale','log')

title('P_{bound} profile')
xlabel('Activator (AU)')
ylabel('P_{bound}')
legend('Original','Fit','Location','northwest')
StandardFigure(gcf,gca)

% Save figures
saveas(gcf,[FigPath,filesep,'P_bound_','parameter_sensitivity_exampleFits' , '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_','parameter_sensitivity_exampleFits' , '.pdf']); 
end