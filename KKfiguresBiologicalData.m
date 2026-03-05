%% Initialize file location
clearvars
% close all

% Extract Current Path
currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

% dataFile{1} = dataFilePath + "Jerome Data\abdomen4.mat";
% dataFile{1} = dataFilePath + "Alex Data\Left Shoulder L123v\ACjoint.mat";
% dataFile{1} = dataFilePath + "Alex Data\Left Shoulder L123v\BicepsTendonXsec.mat";
dataFile{1} = dataFilePath + "Alex Data\Left Shoulder L123v\GlenohumeralJoint.mat";
% dataFile{1} = dataFilePath + "Alex Data\HipBursa.mat";
% dataFile{1} = dataFilePath + "Alex Data\ACJointXsec.mat";

% dataFile{1} = dataFilePath + "KK Data\BiologicalData_2.18.26\LeftThighTop.mat";
% filetype = 2;
filetype = 13;
[p,RFData] = initParams(dataFile,filetype);
p.szAcq = int32(p.szRFframe+1);
%% Process data
% pLarge = computeNewGrid(p,[13,180],[51,350],168*4,300*4);
% pLarge = computeNewGrid(p,[1,p.szX],[1,220],p.szX*2,220*2);
pLarge = computeNewGrid(p,[1,p.szX],[1,300],p.szX*2,300*2);
% pLarge = computeNewGrid(p,[1,p.szX],[1,p.szZ],p.szX*2,p.szZ*2);
% pLarge = p;

M = 7;
tic; images1 = bfmAndProcessFreq(pLarge,RFData,M); toc

M = 19;
tic; images2 = bfmAndProcessFreq(pLarge,RFData,M); toc

%% Plotting
g0 = 0.5;
figBio = plotBiologicalFig_manualPixels(pLarge, images1, images2, g0, 200);
export_fig figBio.png -m4 -transparent


fig = figure('Position',[50 50 612 788]);
plotGammaScaleImage(pLarge.xCoord*1e3,pLarge.zCoord*1e3,images1(1).data,g0)
axis image
export_fig figBioDAS.png -m4 -transparent

%% Helper Functions

function [fig] = plotBiologicalFig_manualPixels(p, images1, images2, g0, input_img_width)

    p.dx = mean(diff(p.xCoord));
    p.dz = mean(diff(p.zCoord));

    % Image aspect ratio in physical units
    aspect_ratio = range(p.zCoord*1e3) / range(p.xCoord*1e3);
    img_height   = input_img_width * aspect_ratio;

    % Grid definition
    num_rows    = 4;
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

    fig = figure('Position',[50 50 fig_width fig_height+30], ...
                 'Color','w');

    % Helper for pixel-position lookup
    tilePosPx = @(n) localTilePositionPixels( ...
        n, num_rows, num_columns, ...
        tile_w, tile_h, gutter_x, gutter_y, fig_height);

    % Calc xOffset
    xOffset = (fig_width-3*tile_w+gutter_x*2)/2;
    
    % Plot the first set of images and label them
    n = 1;
    n = plotSubFigs_manualPixels(fig, tilePosPx, images1, p, n, g0, 30, xOffset); n = n + 1;
    
%     lab2 = annotation(fig, 'textbox', [0,0,.01,.01],...
%         'String', "M=19",'BackgroundColor','white','EdgeColor','none','FontWeight','bold','FontSize',15);
%     pos = tilePosPx(1);
%     posY = pos(2)+30+pos(4);
%     lab2.Units = 'pixels';
%     lab2.Position = [fig_width/2-25,posY,50,30];
    
    % Plot second set of images and label them
    n = plotSubFigs_manualPixels(fig, tilePosPx, images2, p, n, g0, 0, xOffset); n = n + 1;
    
%     lab3 = annotation(fig, 'textbox', [0,0,.01,.01],...
%         'String', "M=48",'BackgroundColor','white','EdgeColor','none','FontWeight','bold','FontSize',15);
%     pos = tilePosPx(8);
%     posY = pos(2);
%     lab3.Units = 'pixels';
%     lab3.Position = [fig_width/2-25,posY,50,30];

end


function [n] = plotSubFigs_manualPixels(fig, tilePosPx, images, p, n, g0, yOffset, xOffset)

    for i = [1, 3:length(images)]

        [~, g] = computeContrastMatch(images(1).data, images(i).data, g0);

        tilePos = tilePosPx(n); tilePos(2) = tilePos(2) + yOffset;
        if (i >= length(images)-1)
            tilePos(1) = tilePos(1) + xOffset;
        end
        
        ax = axes( ...
            'Parent', fig, ...
            'Units',  'pixels', ...
            'Position', tilePos);

        plotGammaScaleImage(p.xCoord*1e3, p.zCoord*1e3, images(i).data, g);
        axis(ax,'image');
        set(ax,'XTick',[],'YTick',[]);
        

        % Label
        x_text     = min(p.xCoord)*1e3 - p.dx*4;
        y_text_top = max(p.zCoord)*1e3 - p.dz*4;

        text(ax, x_text, y_text_top, images(i).name, ...
            'FontSize',12, ...
            'VerticalAlignment','bottom', ...
            'HorizontalAlignment','left', ...
            'Color','black', ...
            'BackgroundColor','white', ...
            'Margin',1);
        
        % Display the top text with black font and white background
        label_str = "M=" + num2str(length(images(i).RXangle(:)));
        text(max(p.xCoord)*1e3 - p.dx*4, max(p.zCoord)*1e3 - p.dz*4, label_str, 'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'black', 'BackgroundColor', 'white', 'Margin', 1, 'Parent', ax)


        n = n + 1;
    end

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
