%% Initialize file location
clearvars
% close all

% Extract Current Path
currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

% dataFile{1} = dataFilePath + "KK Data\TallPhantom_1.5.26\ContrastTargets.mat";
dataFile{1} = dataFilePath + "KK Data\TallPhantom_1.9.26\ContrastTargets_24.mat";
filetype = 2;
[p,RFData] = initParams(dataFile,filetype);
p.szAcq = int32(p.szRFframe+1);
%% Process data
% pLarge = computeNewGrid(p,[71,120],[76,130],50*4,55*4);
pLarge = computeNewGrid(p,[76,140],[71,135],65*4,65*4);
% pLarge = p;

% M = double(pLarge.na);
M = 7;
tic; images = bfmAndProcessFreq(pLarge,RFData,M); toc

M = 19;
tic; images2 = bfmAndProcessFreq(pLarge,RFData,M); toc

%% Plotting

figContrast = plotContrastFig_manualPixels(pLarge, images, images2, 0.5, 200);
% export_fig figContrast.png -m4 -transparent

%% Helper Functions

function [fig] = plotContrastFig_manualPixels(p, images, images2, g0, input_img_width)

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
    n = plotSubFigs_manualPixels(fig, tilePosPx, images,  p, n, g0, 30, xOffset); n = n + 1;
    
    % Plot second set of images and label them
    n = plotSubFigs_manualPixels(fig, tilePosPx, images2, p, n, g0, 0, xOffset); n = n + 1;

end


function [n] = plotSubFigs_manualPixels(fig, tilePosPx, images, p, n, g0, yOffset, xOffset)

    % Define New Region Selector
%     nReg = [141,61;180,170];   % [x1,z1;x2,z2] coordss of opposite corners
%     sReg = [85,34;118,42];     % [x0,rx;z0,rz] center and radii of ellipse
%     iReg = [1,1;p.szX,p.szZ];   % [x1,z1;x2,z2] coords of opposite corners
%     name = "Atten";
    
    nReg = [181,66;250,210];   % [x1,z1;x2,z2] coordss of opposite corners
    sReg = [89,48;138,67];     % [x0,rx;z0,rz] center and radii of ellipse
    iReg = [1,1;p.szX,p.szZ];   % [x1,z1;x2,z2] coords of opposite corners
    name = "Atten";
    
    % Check if area of signal and noise region are roughly equal. if not,
    % throw a warning
    sRegA = round(pi*sReg(3)*sReg(4));
    nRegA = prod(diff(nReg,1)+1);
    
    if ( (nRegA < sRegA*0.9) || (nRegA > sRegA*1.1) )
        warning("Number of pixels in signal Region and noise region are not in the same ballpark! " + ...
            "nRegPixels = " + num2str(nRegA) + ". sRegPixels = " + num2str(sRegA) + ".");
    end
    
    Rnew = RegionSelector(p,nReg,sReg,iReg,name);
    

    imgPmatch = repmat(struct("data",[],"name",[]),1,length(images));
    for i = 1:length(images)
        % Extract Matched images on a per frame basis
        [imgPmatch(i).data,g(n,i)] = computeContrastMatch(images(1).data, images(i).data, g0);
    end
    
    % Compute Metrics for current frame
    [~, ~, GCNR] = computeContMetrics(Rnew,imgPmatch);
    


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
        
        if (i==1)

            hold(ax,'on');

            % ---- Red rectangle (pixel coords → physical coords) ----
            xPix = [nReg(1,1), nReg(2,1)];
            zPix = [nReg(1,2), nReg(2,2)];

            xVals = p.xCoord(xPix) * 1e3;
            zVals = p.zCoord(zPix) * 1e3;

            rectPos = [ ...
                min(xVals), ...
                min(zVals), ...
                abs(diff(xVals)), ...
                abs(diff(zVals)) ];

            rectangle(ax, ...
                'Position', rectPos, ...
                'EdgeColor','r', ...
                'LineWidth',1.5);

            % ---- Ellipse (pixel center + radii → physical coords) ----
            x0  = p.xCoord(sReg(1,1)) * 1e3;
            z0  = p.zCoord(sReg(2,1)) * 1e3;

            rx  = abs(p.xCoord(sReg(1,1)+sReg(1,2)) ...
                    - p.xCoord(sReg(1,1))) * 1e3;

            rz  = abs(p.zCoord(sReg(2,1)+sReg(2,2)) ...
                    - p.zCoord(sReg(2,1))) * 1e3;

            ellipsePos = [ ...
                x0 - rx, ...
                z0 - rz, ...
                2*rx, ...
                2*rz ];

            rectangle(ax, ...
                'Position', ellipsePos, ...
                'Curvature',[1 1], ...
                'EdgeColor','b', ...
                'LineWidth',1.5);

            hold(ax,'off');

        end

        % Label
        x_text     = min(p.xCoord)*1e3 - p.dx*4;
        y_text_top = min(p.zCoord)*1e3;

        text(ax, x_text, y_text_top, images(i).name, ...
            'FontSize',12, ...
            'VerticalAlignment','top', ...
            'HorizontalAlignment','left', ...
            'Color','black', ...
            'BackgroundColor','white', ...
            'Margin',1);
        
        % Display the top text with black font and white background
        label_str = "M=" + num2str(length(images(i).RXangle(:)));
        text(max(p.xCoord)*1e3 - p.dx*4, min(p.zCoord)*1e3, label_str, 'FontSize', 12, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', 'black', 'BackgroundColor', 'white', 'Margin', 1, 'Parent', ax)

        
        % Display the bottom text with black font and white background
        label_str = sprintf('%.3f', GCNR(i));
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


function [C, CNR, GCNR] = computeContMetrics(R,images)

    nImg = length(images);
    signalMap = R.getsMap(1);
    noiseMap = R.getnMap(1);
    imgMap = R.getiMap(1);
    C = zeros(nImg,1); CNR = C; GCNR = C;

    for n = 1:nImg
        img = images(n).data;

        mu_i=mean(img(signalMap)); mu_o=mean(img(noiseMap));
        v_i=var(img(signalMap)); v_o=var(img(noiseMap));

        C(n)=10*log10(mu_i./mu_o);
        CNR(n)=10*log10(abs(mu_i-mu_o)/sqrt(v_i+v_o));

        % Select subset of image to generate histogram for. Compute w.r.t.
        % magnitudes
        imgReg = img(imgMap);
        x=linspace(min(imgReg(:)),max(imgReg(:)),25);


        % [~,edges] = histcounts(img(signalMap), 'Normalization','pdf');
        % pdf_i = edges(2:end) - (edges(2)-edges(1))/2;
        % 
        % [~,edges] = histcounts(img(noiseMap), 'Normalization','pdf');
        % pdf_o = edges(2:end) - (edges(2)-edges(1))/2;

        [pdf_i]=hist(img(signalMap),x);
        [pdf_o]=hist(img(noiseMap),x);

        OVL=sum(min([pdf_i./sum(pdf_i); pdf_o./sum(pdf_o)]));
        GCNR(n)= 1 - OVL;

    end

end
