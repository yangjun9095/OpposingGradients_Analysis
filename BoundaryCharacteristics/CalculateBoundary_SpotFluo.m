function [Xhalf,Width] = CalculateBoundary_SpotFluo(Pattern,Time,nc12,nc13,nc14)
% Description : This function plots the pattern (mRNA/protein) over AP, at
% different time points, and shows boundary position, width, etc. 

% INPUT : Pattern of Spot fluorescence/mRNA/protein over AP, over time (Time x APbins)

% OUTPUT : Xhalf (boundary position), and Width (boundary width), which has
% dimension of (Time x 1).


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
        SmoothPattern(j,i) = nanmean(Pattern(j,i-2:i+2));
    end
end

Pattern = SmoothPattern;

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
    
    title({'MS2 Spot','Fluorescence'})
    xlabel('AP')
    ylabel('Fluorescence(AU)')
    xlim([0.15 0.6])
    ylim([0 1000])
    drawnow
    hold off
    
    %pchip fitting
    dxx=0.001;
    xx=0.2:dxx:0.8;
    yy=pchip(APbin,Pattern(tpoint,:),xx);
    
    %Find the x position where it has half of the Maximum intensity &
    %negative slope (find the most anterior one)
    if (sum(yy<=0))==length(yy)
        Xhalf(tpoint)=nan;
        HalfMaxY(tpoint)=nan;
        Index(tpoint)=nan;
    else
        i=11;
        while ~(sum(yy<=0)==length(yy))
            if yy(i)<=0.5*MaxIntensity(tpoint)&&...
                yy(i-1)>0.5*MaxIntensity(tpoint)&&...
                (yy(i+10)-yy(i-10))/(xx(i+10)-xx(i-10)) <0
                Xhalf(tpoint)=xx(i);
                HalfMaxY(tpoint)=yy(i);
                Index(tpoint)=i;
                break
            else 
                i=i+1;
            end
        end
    end
        
            
    
%     for i=11:length(xx)-11
%         if (yy(i)<=0.5*MaxIntensity(tpoint)) && ~(sum(yy<=0)==length(yy))
%             if yy(i-1)>0.5*MaxIntensity(tpoint) && (yy(i+10)-yy(i-10))/(xx(i+10)-xx(i-10)) <0
%                 Xhalf(tpoint)=xx(i);
%                 HalfMaxY(tpoint)=yy(i);
%                 Index(tpoint)=i;
%             end
%         elseif (yy(i)<=0.5*MaxIntensity(tpoint)) && (sum(yy<=0))==length(yy) % In case yy pattern is all flat zero.
%                 Xhalf(tpoint)=nan;
%                 HalfMaxY(tpoint)=nan;
%                 Index(tpoint)=nan;
%         end
%     end
    

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
    if ~(sum(yy<=0)==length(yy))
        if tpoint<nc13 % for earlier cycles, the pattern 
            for i=26:length(yy)-26
                grad(i)=(yy(i+25)-yy(i-25)) / (50*dxx);
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

        Width(tpoint) = (Right-Left)*dX;

        
    else
        Width(tpoint) = nan;
        Y=zeros(size(X));
        Maximum = 0;
        Minimum = 0;
        Left = 0;
        Right = 0;
    end

    
    subplot(1,4,2)
    plot(APbin,Pattern(tpoint,:),'o','color',Color(tpoint-iStart+1,:))
    hold on
    plot(xx,yy,'color',Color(tpoint-iStart+1,:))
    xlim([0.2 0.6])
    %ylim([0 30000])
    title({'Boundary','Calculation'})
    xlabel('AP')
    ylabel('Spot fluorescence(AU)')
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
    ylabel('Spot fluo Boundary')
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
    ylabel('Spot Fluo Width')
    hold off
    %pause
    %saveas(gcf,['E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Movies-OpposingGradients\BoundaryFeature-r3\3APbinsAveraged\BoundaryFeature','r3 @ ','Time = ',num2str(floor(Time(tpoint))),'min from nc12'],'tiff')

end

%v = VideoWriter('mRNA.avi');



end