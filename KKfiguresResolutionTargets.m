%% KKfiguresResolutionTargets - Generate resolution target analysis figures
%
% This script generates publication-quality figures comparing DAS and KK
% beamforming on resolution phantom targets. Includes zoom regions and
% lateral line profiles for PSF width assessment.
%
% User paths to modify:
%   - dataFilePath: path to ultrasound dataset directory
%
% Required data: ResolutionTargets .mat files
%
% Required functions: initParams, bfmAndProcessFreq, computeNewGrid,
%   computeContrastMatch, plotGammaScaleImage
%
% Outputs: Full-view figures, zoom figures, and lateral intensity profiles

%% Initialize file location
clearvars

% Extract Current Path
currentDir = matlab.desktop.editor.getActiveFilename;
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

dataFile{1} = dataFilePath + "KK Data\TallPhantom_12.17.25\ResolutionTargets_48.mat";
filetype = 0;
[p,RFData] = initParams(dataFile,filetype);
p.szAcq = int32(p.szRFframe+1);
%% Process data
pLarge = computeNewGrid(p,[51,140],[141,210],90*4,70*4);

M = double(pLarge.na);
tic; images = bfmAndProcessFreq(pLarge,RFData,M); toc
%% Plotting
zC = [191,235]; xC = [31,90];

%% Plotting full data
figure('Units','pixels','Position',[357,68,1337,909])
tiledlayout('flow')

for i = 1:length(images)
    nexttile
    plotGammaScaleImage(images(i).data,0.5)
    title(images(i).name);
end
%% Helper Functions

function [fig] = plotResFig(p, images, g0, input_img_width, xC, zC)

    p.dx = mean(diff(p.xCoord));
    p.dz = mean(diff(p.zCoord));

    % Calculate the aspect ratio of the first image
    aspect_ratio = range(p.zCoord*1e3) / range(p.xCoord*1e3);
    
    % Calculate the height of each image based on the aspect ratio
    img_height = input_img_width * aspect_ratio;

    % Define the number of rows and columns
    num_rows = 6;
    num_columns = 4;

    % Calculate the total figure width and height based on the number of rows and columns
    scale = 1;
    fig_width = input_img_width * num_columns * scale;
    fig_height = img_height * num_rows * scale;

    % Create the figure with the calculated dimensions
    fig = figure('Position', [50, 50, fig_width, fig_height+5]);  % [left, bottom, width, height]

    t = tiledlayout(num_rows,num_columns,'TileSpacing','tight','Padding','none');
    n = 1;
    
    for i = [1,3:length(images)]
        
        [~,g] = computeContrastMatch(images(1).data,images(i).data,g0);
        
        nexttile(n,[1,1])
        plotGammaScaleImage(p.xCoord*1e3,p.zCoord*1e3,images(i).data,g)

        axis image
        set(gca, 'xtick', [], 'ytick', []);

        % Plot Zoom Box on DAS image
        if (i==1)
            % Convert regCoords to millimeters
            xCoords_mm = p.xCoord(xC)*1e3;
            zCoords_mm = p.zCoord(zC)*1e3;

            % Define the rectangle position and size
            rectX = xCoords_mm(1);
            rectY = zCoords_mm(1);
            rectWidth = xCoords_mm(2) - xCoords_mm(1);
            rectHeight = zCoords_mm(2) - zCoords_mm(1);

            % Plot the rectangle
            hold on;
            rectangle('Position', [rectX, rectY, rectWidth, rectHeight], 'EdgeColor', 'r', 'LineWidth', 2);
            hold off;
        end
        
        % Add the subfigure label
        label_str = images(i).name;
        
        % Text position
        x_text = min(p.xCoord)*1e3 - p.dx*4;
        y_text_top = min(p.zCoord)*1e3;
        
        % Display the top text with black font and white background
        text(x_text, y_text_top, label_str, 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'black', 'BackgroundColor', 'white', 'Margin', 1, 'Parent', gca)
        n = n+1;
    end
    
    img2 = images;
    for i = 1:length(images)
        img2(i).data = images(i).data(zC(1):zC(2),xC(1):xC(2));
    end
    pZ = computeNewGrid(p,[xC(1),xC(2)],[zC(1),zC(2)]);
    
    n = n+1;
    for i = [1,3:length(images)]
        
        [~,g] = computeContrastMatch(img2(1).data,img2(i).data,g0);
        
        nexttile(n,[1,1])
        plotGammaScaleImage(pZ.xCoord*1e3,pZ.zCoord*1e3,img2(i).data,g)

        axis image
        set(gca, 'xtick', [], 'ytick', []);
        
        % Plot a red dashed line at this y-coordinate
        if (i == 1)
            hold on;
            plot(pZ.xCoord*1e3, repmat(pZ.zCoord(23)*1e3, size(pZ.xCoord*1e3)), 'r--', 'LineWidth', 2);
            hold off;
        end

        % Add the subfigure label
        label_str = images(i).name + " zoom";
        
        % Text position
        x_text = min(pZ.xCoord)*1e3 - pZ.dx*4;
        y_text_top = min(pZ.zCoord)*1e3;
        
        % Display the top text with black font and white background
        text(x_text, y_text_top, label_str, 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'black', 'BackgroundColor', 'white', 'Margin', 1, 'Parent', gca)
        n = n+1;
    end

    lPf = zeros(p.szX,4);
    legLabels = cell(4,1);
    k = 1;
    for i = [1,3:5]
        [~,g] = computeContrastMatch(images(1).data,images(i).data,g0);
        slice = mean(abs(images(i).data(213:214,:)).^g,1);
        lPf(:,k) = slice./max(slice(:));
        legLabels{k} = images(i).name; k = k+1;
    end
    
    lcolors = {"b","r","m","g"};

    nexttile(n+1,[2,num_columns])
    hold on
    for i = 1:k-1
        pline = plot(p.xCoord*1e3,lPf(:,i));
        pline.Color = lcolors{i};
    end
    legend(legLabels{:})
    set(gca,'xlim',[p.xCoord(1),p.xCoord(end)]*1e3,'Box','off')
    xlabel('Lateral Position (mm)')
    ylabel('Normalized Intensity')
end

function [fig] = plotResFig_manualPosition(p, images, g0, input_img_width, xC, zC)

    p.dx = mean(diff(p.xCoord));
    p.dz = mean(diff(p.zCoord));

    % Calculate the aspect ratio of the first image
    aspect_ratio = range(p.zCoord*1e3) / range(p.xCoord*1e3);
    
    % Calculate the height of each image based on the aspect ratio
    img_height = input_img_width * aspect_ratio;

    % Define the number of rows and columns
    num_rows = 6;
    num_columns = 4;
    
    % Pixel gutters (set to 0 for fully flush)
    gutter_x = 2;
    gutter_y = 2;

    % Tile size in pixels
    tile_w = input_img_width;
    tile_h = img_height;

    % Figure size
    fig_width  = num_columns * tile_w + (num_columns-1)*gutter_x;
    fig_height = num_rows    * tile_h + (num_rows-1)*gutter_y;

    % Create the figure with the calculated dimensions
    fig = figure('Position', [50, 50, fig_width, fig_height+30]);  % [left, bottom, width, height]

    tilePosPx = @(num) localTilePositionPixels( ...
        num, num_rows, num_columns, ...
        tile_w, tile_h, gutter_x, gutter_y, fig_height);
    
    % Calc xOffset
    xOffset = (fig_width-3*tile_w+gutter_x*2)/2;
    
    n = 1;
    
    for i = [1,3:length(images)]
        
        [~,g] = computeContrastMatch(images(1).data,images(i).data,g0);

        tilePos = tilePosPx(n);
        tilePos(2) = tilePos(2) + 30;
        if (i >= length(images)-1)
            tilePos(1) = tilePos(1) + xOffset;
        end
        ax = axes( ...
            'Parent', fig, ...
            'Units',  'pixels', ...
            'Position', tilePos);
        
        
        plotGammaScaleImage(ax,p.xCoord*1e3,p.zCoord*1e3,images(i).data,g)

        axis image
        set(gca, 'xtick', [], 'ytick', []);

        % Plot Zoom Box on DAS image
        if (i==1)
            % Convert regCoords to millimeters
            xCoords_mm = p.xCoord(xC)*1e3;
            zCoords_mm = p.zCoord(zC)*1e3;

            % Define the rectangle position and size
            rectX = xCoords_mm(1);
            rectY = zCoords_mm(1);
            rectWidth = xCoords_mm(2) - xCoords_mm(1);
            rectHeight = zCoords_mm(2) - zCoords_mm(1);

            % Plot the rectangle
            hold on;
            rectangle('Position', [rectX, rectY, rectWidth, rectHeight], 'EdgeColor', 'r', 'LineWidth', 2);
            hold off;
        end
        
        % Display the top text with black font and white background
        label_str = "M=" + num2str(length(images(i).RXangle(:)));
        text(max(p.xCoord)*1e3 - p.dx*4, min(p.zCoord)*1e3, label_str, 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', 'black', 'BackgroundColor', 'white', 'Margin', 1, 'Parent', ax)
        
        % Add the subfigure label
        label_str = images(i).name;
        
        % Text position
        x_text = min(p.xCoord)*1e3 - p.dx*4;
        y_text_top = min(p.zCoord)*1e3;
        
        % Display the top text with black font and white background
        text(x_text, y_text_top, label_str, 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'black', 'BackgroundColor', 'white', 'Margin', 1, 'Parent', gca)
        n = n+1;
    end
    
    img2 = images;
    for i = 1:length(images)
        img2(i).data = images(i).data(zC(1):zC(2),xC(1):xC(2));
    end
    pZ = computeNewGrid(p,[xC(1),xC(2)],[zC(1),zC(2)]);
    
    n = n+1;
    for i = [1,3:length(images)]
        
        [~,g] = computeContrastMatch(img2(1).data,img2(i).data,g0);
        
        tilePos = tilePosPx(n);
        tilePos(2) = tilePos(2) + 20;
        if (i >= length(images)-1)
            tilePos(1) = tilePos(1) + xOffset;
        end
        ax = axes( ...
            'Parent', fig, ...
            'Units',  'pixels', ...
            'Position', tilePos);
        
        plotGammaScaleImage(ax, pZ.xCoord*1e3,pZ.zCoord*1e3,img2(i).data,g)

        axis image
        set(ax, 'xtick', [], 'ytick', []);
        
        % Plot a red dashed line at this y-coordinate
        if (i == 1)
            hold on;
            plot(pZ.xCoord*1e3, repmat(pZ.zCoord(23)*1e3, size(pZ.xCoord*1e3)), 'r--', 'LineWidth', 2);
            hold off;
        end

        % Add the subfigure label
        label_str = images(i).name + " zoom";
        
        % Text position
        x_text = min(pZ.xCoord)*1e3 - pZ.dx*4;
        y_text_top = min(pZ.zCoord)*1e3;
        
        % Display the top text with black font and white background
        text(x_text, y_text_top, label_str, 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Color', 'black', 'BackgroundColor', 'white', 'Margin', 1, 'Parent', gca)
        n = n+1;
    end

    lPf = zeros(p.szX,4);
    legLabels = cell(4,1);
    k = 1;
    for i = [1,3:5]
        [~,g] = computeContrastMatch(images(1).data,images(i).data,g0);
        slice = mean(abs(images(i).data(213:214,:)).^g,1);
        lPf(:,k) = slice./max(slice(:));
        legLabels{k} = images(i).name; k = k+1;
    end
    
    lcolors = {"b","r","m","g"};

    tilePos = tilePosPx(n+1+num_columns);
    tilePos(3) = tilePos(3)*num_columns - 40;
    tilePos(1) = tilePos(1) + 40;
    tilePos(4) = tilePos(4)*2 - 20;
    tilePos(2) = tilePos(2) + 40;
        
    ax = axes( ...
        'Parent', fig, ...
        'Units',  'pixels', ...
        'Position', tilePos);
    hold(ax,'on')
    for i = 1:k-1
        pline = plot(ax,p.xCoord*1e3,lPf(:,i));
        pline.Color = lcolors{i};
    end
    legend(ax,legLabels{:})
    set(ax,'xlim',[p.xCoord(1),p.xCoord(end)]*1e3,'Box','off')
    xlabel(ax,'Lateral Position (mm)')
    ylabel(ax,'Normalized Intensity')
end

function pos = localTilePositionPixels(n, num_rows, num_cols, ...
                                       tile_w, tile_h, gutter_x, gutter_y, fig_h)
    % Row-major indexing (like nexttile)
    row = ceil(n / num_cols);           % 1 = top row
    col = n - (row-1)*num_cols;

    % Left coordinate
    x = (col-1) * (tile_w + gutter_x);

    % Convert top-based row index to MATLAB bottom-based pixels
    y_top = (row-1) * (tile_h + gutter_y);
    y     = fig_h - y_top - tile_h;

    pos = [x, y, tile_w, tile_h];
end

