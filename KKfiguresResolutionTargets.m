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
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"KKBeamforming"),1)},"Datasets");

dataFile{1} = fullfile(dataFilePath, "ResolutionTargets_48.mat");
filetype = 0;
[p,RFData] = initParams(dataFile,filetype);
p.szAcq = int32(p.szRFframe+1);
%% Process data
pLarge = computeNewGrid(p,[51,140],[141,210],90*4,70*4);

M = double(pLarge.na);
tic; images = bfmAndProcessFreq(pLarge,RFData,M); toc
%% Plotting
zC = [191,235]; xC = [31,90];

figRes = plotResFig(pLarge,images,0.5,200,xC,zC);

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

    lcolors = {"k","b","r","g"};
    xRng = 26:96;
    n = n+1;
    plotLatResSubplot(images,g0,p,lcolors,xRng,n,false,true)

    xRng = 100:130;
    n = n+1;
    plotLatResSubplot(images,g0,p,lcolors,xRng,n,false,false)

    xRng = 154:184;
    n = n+1;
    plotLatResSubplot(images,g0,p,lcolors,xRng,n,false,false)

    xRng = 226:256;
    n = n+1;
    plotLatResSubplot(images,g0,p,lcolors,xRng,n,true,false)
    
end


function plotLatResSubplot(images,g0,p,lcolors,xRng,n,lFlag,yFlag)


    lPf = zeros(length(xRng),4);
    legLabels = cell(4,1);
    k = 1;
    for i = [1,3:5]
        [~,g] = computeContrastMatch(images(1).data,images(i).data,g0);
        slice = mean(abs(images(i).data(213:214,xRng)).^g,1);
        lPf(:,k) = slice./max(slice(:));
        legLabels{k} = images(i).name; k = k+1;
    end


    nexttile(n,[2,1])
    hold on
    for i = 1:k-1
        pline = plot(p.xCoord(xRng)*1e3,lPf(:,i));
        pline.Color = lcolors{i};
        if i == 1
            pline.LineStyle = "--";
        end
    end
    if (lFlag)
        legend(legLabels{:})
    end
    set(gca,'xlim',[p.xCoord(xRng(1)),p.xCoord(xRng(end))]*1e3,'Box','off')
    xlabel('Lateral Position (mm)')

    if (yFlag)
        ylabel('Normalized Intensity')
    else
        set(gca,'ytick',[]);
    end



end

