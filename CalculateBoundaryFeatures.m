function [Xhalf,Width] = CalculateBoundaryFeatures(Pattern,Time,nc12,nc13,nc14)
% Description : This function plots the pattern (mRNA/protein) over AP, at
% different time points, and shows boundary position, width, etc. 

% INPUT : Pattern of mRNA/protein over AP, over time (Time x APbins)

% OUTPUT : 

%% Load the dataset
%Data = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\r1_FromNC12.mat')
%Pattern = Data.AccumulatedmRNA_FractionON;
%Time = Data.ElapsedTime;
% Define fields for future use
%nc12 = Data.nc12;
%nc13 = Data.nc13;
%nc14 = Data.nc14;

%% Smoothening the pattern
% Smoothen the Pattern over 3 AP windows, neighboring AP bins (7.5%)
SmoothPattern = nan(length(Pattern),41);
for j=1:length(Pattern) %time length
    for i=8:33 % 20%~80% AP axis
        SmoothPattern(j,i) = nanmean(Pattern(j,i-1:i+1));
    end
end

Pattern = SmoothPattern;
%% Analyze the boundary position / width for AccumulatedmRNA
% %In this dataset, I need to consider only 0.2~0.6 since it's unlikely that boundary exist outside this range.
% %Input : Pattern (time x APbin), assuming that we've already smoothened
% %the pattern.
% %For specific time point, Pattern(tpoint,:) -> mRNA at that time point.
% APbin=[0:0.025:1];
% tpoint=14;
% 
% % plot(0:0.025:1,Pattern(tpoint,:))
% % xlim([0.25 0.6])
% hold on
% %Spline fitting
% dxx=0.001;
% xx=0.2:dxx:0.6;
% yy=spline(APbin,Pattern(tpoint,:),xx);
% yyy=pchip(APbin,Pattern(tpoint,:),xx);
% yyyy=interp1(APbin,Pattern(tpoint,:),xx);
% %figure(2)
% hold on
% plot(APbin,Pattern(tpoint,:),'o',xx,yy,'r')
% plot(xx,yyy,'b')    %APbin,Pattern(tpoint,:),'o',
% plot(xx,yyyy,'g')  %APbin,Pattern(tpoint,:),'o'
% xlim([0.2 0.6])
% 
% legend('data','spline','pchip','interp')

% %% Find the x position where it has half of the Maximum intensity
% % Making Nans/inf to zeros, for maximum
% Pattern(isnan(Pattern)) = 0;
% Pattern(Pattern==inf) = 0;
% 
% for t=1:length(Pattern)
%     MaxIntensity(t)=max(Pattern(t,:));
% end
% 
% for i=1:length(xx)
%     if yy(i)<0.5*MaxIntensity(tpoint)
%         if yy(i-1)>=0.5*MaxIntensity(tpoint)
%             Xhalf=xx(i);
%             HalfMaxY=yy(i);
%             Index(tpoint)=i;
%         end    
% %     else 
% %         error('sth about pattern is wrong')
%     end
% end
% %Xhalf
% %need to think about calculating the width.
% line([Xhalf Xhalf],[0 HalfMaxY])
% %ylim([0 19000])
% title('Fitting Curve to the Data')
% xlabel('AP')
% ylabel('mRNA(AU)')
% set(gca,'fontsize',30)
% hold off
%% Plot for AccumulatedmRNA and boundary position/width for all time points

iStart=2; %Start time
iEnd=length(Pattern); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
    
% Find the maximum value of pattern for each time point

% Making Nans/inf to zeros, for maximum
Pattern(isnan(Pattern)) = 0;
Pattern(Pattern==inf) = 0;

for t=1:length(Pattern)
    MaxIntensity(t)=max(Pattern(t,:));
end

% Calculate the Boundary position/width
for tpoint=2:length(Pattern)
    
    APbin=[0:0.025:1];

    % plot(0:0.025:1,Pattern(tpoint,:))
    % xlim([0.25 0.6])
    subplot(1,4,1)
    hold on
    plot(APbin,Pattern(tpoint,:),'color',Color(tpoint-iStart+1,:))
    
    title({'mRNA','Accumulation'})
    xlabel('AP')
    ylabel('mRNA Accumulation(AU)')
    xlim([0.15 0.6])
    ylim([0 30000])
    drawnow
    hold off
    
    %pchip fitting
    dxx=0.001;
    xx=0.2:dxx:0.6;
    yy=pchip(APbin,Pattern(tpoint,:),xx);
    
    %Find the x position where it has half of the Maximum intensity &
    %negative slope
    for i=26:length(xx)-26
        if (yy(i)<=0.5*MaxIntensity(tpoint))
            if yy(i-1)>0.5*MaxIntensity(tpoint) && (yy(i+25)-yy(i-25))/(xx(i+25)-xx(i-25)) <0
                Xhalf(tpoint)=xx(i);
                HalfMaxY(tpoint)=yy(i);
                Index(tpoint)=i;
            end
        elseif (yy(i)<=0.5*MaxIntensity(tpoint)) && (sum(yy==0))==length(yy) % In case yy pattern is all flat zero.
        end
    end
    
    %Get the Width using xx,yy pchip fit
    %First, draw a tangent line at the midpoint.
    %Mid point
    MidX(tpoint)=Xhalf(tpoint);
    MidY(tpoint)=HalfMaxY(tpoint);
    
    %tangent line
    dX=0.001;
    X=0:dX:1;
    % Get the slope (grad) from 2.5% neighboring APbin (toward anterior,
    % this is for convenience, needs justification later)
    if tpoint<nc13 % for earlier cycles, the pattern 
        for i=11:length(yy)-11
            grad(i)=(yy(i+10)-yy(i-10)) / (20*dxx);
        end
    else
        for i=26:length(yy)-26
            grad(i)=(yy(i+25)-yy(i-25)) / (50*dxx);
        end
    end
        
    %grad=diff(yy)./diff(xx);
    %grad(Index) : slope(gradient) at the midpoint.
    MidIndex=Index(tpoint);
    Y=grad(MidIndex)*(X-MidX(tpoint))+MidY(tpoint);
    
    Maximum=max(yy);
    Minimum=min(yy);
    
    for i=1:length(X)-1
        if Y(i)>Maximum
            if Y(i+1)<Maximum
                Left=i;
            %else
                %print('sth went wrong')
            end
        end
    end
    
    for j=2:length(X)
        if Y(j)<=Minimum
            if Y(j-1)>=Minimum
                Right=j;
            %else
                %print('sth went wrong')
            end
        end
    end
    
    % Get the Width only when the pattern is non-zero
    if (yy==0)
        Right = 0;
        Left = 0;
        Width(tpoint) = 0;
    else
        Width(tpoint) = (Right-Left)*dX;
    end
%     tpoint
%     Xhalf(tpoint)
    
    subplot(1,4,2)
    plot(APbin,Pattern(tpoint,:),'o','color',Color(tpoint-iStart+1,:))
    hold on
    plot(xx,yy,'color',Color(tpoint-iStart+1,:))
    xlim([0.2 0.6])
    %ylim([0 30000])
    title({'Boundary','Calculation'})
    xlabel('AP')
    ylabel('mRNA Accumulation(AU)')
    %hold off
    
    %need to think about calculating the width.
    line([Xhalf(tpoint) Xhalf(tpoint)],[0 HalfMaxY(tpoint)])
    plot(X,Y,'k')
    bar(Maximum)
    bar(Minimum)
    line([Left*dX Right*dX],[Maximum+1 Maximum+1])
    
    hold off
    %drawnow
    %pause
    
    
    subplot(1,4,3)
    hold on
    plot(Time(tpoint),Xhalf(tpoint),'o','color',Color(tpoint-iStart+1,:))
    line([Time(nc13) Time(nc13)],[0 0.5])
    line([Time(nc14) Time(nc14)],[0 0.5])
    xlim([0 max(Time)+1])
    ylim([0.2 0.5])
    title({'Boundary Position','along Time'})
    xlabel('Time (min)')
    ylabel('mRNA Boundary')
    hold off
    %pause(0.5)
    
    subplot(1,4,4)
    hold on
    plot(Time(tpoint),Width(tpoint),'o','color',Color(tpoint-iStart+1,:))
    line([Time(nc13) Time(nc13)],[0 1])
    line([Time(nc14) Time(nc14)],[0 1])
    xlim([0 max(Time)+1])
    ylim([0 1])
    
    title({'Boundary Width','along Time'})
    xlabel('Time (min)')
    ylabel('mRNA Boundary Width')
    hold off
    %pause
    %saveas(gcf,['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Movies-OpposingGradients\BoundaryFeature-r3\3APbinsAveraged\BoundaryFeature','r3 @ ','Time = ',num2str(floor(Time(tpoint))),'min from nc12'],'tiff')

end

%v = VideoWriter('mRNA.avi');



end