function standardizeFigure(ax, legend, varargin)

try
    colorDict = struct();
    colorDict.brown = [179.155,142]/255;
    colorDict.red = [213,108,85]/255;
    colorDict.yellow = [234,194,100]/255;
    colorDict.cyan = [108,188,233]/255;
    colorDict.magenta = [208,109,171]/255;
    colorDict.lightBlue = [115,142,193]/255;
    colorDict.blue = [115,143,193]/255;
    colorDict.purple = [171,133,172]/255;
    colorDict.green = [122,169,116]/255;
    colorDict.darkgreen = [126,157,144]/255;
    
    colorDictFields = fields(colorDict);
    
    color(1,:) = [0 0 0];

    axesLineWidth = 5;
    fig = gcf;
    legend = findobj(fig, 'Type', 'Legend');
    dataObj = get(ax, 'Children');
    dataType = get(dataObj, 'Type');
    if ~iscell(dataType)
        dataType = {dataType};
    end
    legendSize = 20;
    fontSize = 20;
    
    for i = 1:length(varargin)
       if strcmpi(varargin{i}, 'axeslinewidth')
            axesLineWidth = varargin{i+1};            
        elseif strcmpi(varargin{i}, 'red')
            color(i,:) = [213,108,85]/255;
        elseif strcmpi(varargin{i}, 'yellow')
            color(i,:) = [234,194,100]/255;
        elseif strcmpi(varargin{i}, 'cyan')
            color(i,:) = [108,188,233]/255;
        elseif strcmpi(varargin{i}, 'magenta')
            color(i,:) = [208,109,171]/255;
        elseif strcmpi(varargin{i}, 'lightBlue')
            color(i,:) = [115,142,193]/255;
        elseif strcmpi(varargin{i}, 'legendFontSize')
            legendSize = varargin{i+1};
        elseif strcmpi(varargin{i}, 'fontSize')
            fontSize = varargin{i+1};
        end
    end
    
    if ~isempty(legend)
        legend.FontSize = legendSize;
        legend.Box = 'off';
    end
 
    for i = 1:length(dataObj)
        if strcmpi(dataType{i}, 'scatter')
            dataObj(i).Marker = '.';
            dataObj(i).SizeData = 250;
            %Change color to physical biology colors as long as the number
            %of colors needed is less than 5.
            if i <= length(colorDictFields)
%                 dataObj(i).Color = colorDict.(colorDictFields{i});
                dataObj(i).MarkerFaceColor = colorDict.(colorDictFields{i});
                dataObj(i).MarkerEdgeColor = colorDict.(colorDictFields{i});
            end
        elseif strcmpi(dataType{i}, 'bar') || strcmpi(dataType{i}, 'histogram')
%             dataObj(i).LineStyle = 'none';
            if i <= length(colorDictFields)
                dataObj(i).FaceColor = colorDict.(colorDictFields{i});
            end
        elseif strcmpi(dataType{i}, 'line') || strcmpi(dataType{i}, 'errorbar')
            dataObj(i).LineWidth = 5;
            dataObj(i).Marker = '.';
            dataObj(i).MarkerSize = 30;
            %Change color to physical biology colors as long as the number
            %of colors needed is less than 5.
            if i <= length(colorDictFields)
                dataObj(i).Color = colorDict.(colorDictFields{i});
                dataObj(i).MarkerFaceColor = colorDict.(colorDictFields{i});
                dataObj(i).MarkerEdgeColor = colorDict.(colorDictFields{i});
            end
            if strcmpi(dataType{i}, 'errorbar')
                %insert errorbar specific things here.
            end
        end
    end
    
    set(ax, 'TickLength',[0 0],...
        'FontSize', fontSize, 'FontName', 'Myriad Pro', 'FontWeight', 'bold');
    ax.TickDir = 'out';
    ax.LineWidth = axesLineWidth;
    faceColor = [255,251,206]/255; %yellow axis face.
    ax.Color = faceColor;
    fig.Color = [255,255,255]/255; %white figure background
end
end