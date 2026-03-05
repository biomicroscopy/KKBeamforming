clearvars
close all


currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
nL = length(currentDir);

ipath = ['-I' fullfile(currentDir{1:nL-2},'Verasonics Integrated C Code','lib','eigen-3.4.0')]; % mfilename also includes the filename in the directory, so we need to remove that
%% Compile
mingwFlags = {'CXXFLAGS="$CXXFLAGS -march=native -std=c++14 -fno-math-errno -ffast-math -fopenmp -DNDEBUG -w -Wno-error"',...
            'LDFLAGS="$LDFLAGS -fopenmp"','CXXOPTIMFLAGS="-O3"'};
        
% tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'DASLUT.cpp'); toc
% tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'DASLUTV2.cpp'); toc
% tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'DASLUTV3.cpp'); toc
tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'DASLUTV4.cpp'); toc
% tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'DASLUTV5.cpp'); toc

% tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'RXLUTTest.cpp'); toc

%% Load Data
currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

dataFile{1} = dataFilePath + "KK Data\TallPhantom_12.17.25\ResolutionTargets_48.mat";
filetype = 2;

% Load and initialize data and parameter struct
[p,RFData] = initParams(dataFile,filetype);
p.szAcq = int32(p.szRFframe+1);

% Reorganize data into [samples x channels x angles] and hilbert transform
cRF = computeCRF(double(RFData),p); 

%% Initialize
pL = initParamsLUT(p);
pL2 = initParamsLUTV2(p);

pL = convertParamsToSingle(pL);
pL = convertParamsToInteger(pL);
pL2 = convertParamsToSingle(pL2);
pL2 = convertParamsToInteger(pL2);

beamform = reconraw.DASBModeOffline(p);

cRF2 = UnscrambleRF(pL,cRF);
cRF2 = single(reshape(cRF2,[p.szAcq,p.numEl*p.na]));
cRF3 = single(cRF);
%% Execute


tic;
Recon = DASLUT(pL,cRF2,pL.RXDelay,pL.TXDelay);
toc

tic;
Recon2 = DASLUTV2(pL2,cRF2,pL2.RXDelayX,pL2.RXDelayZ,pL2.TXDelayX,pL2.TXDelayZ);
toc

tic;
Recon3 = DASLUTV3(pL2,cRF2,pL2.RXDelay,pL2.TXDelayX,pL2.TXDelayZ);
toc

TXDelayX = pL2.TXDelayX.';
TXDelayZ = pL2.TXDelayZ.';

tic;
Recon4 = DASLUTV4(pL2,cRF3,pL2.RXDelay,TXDelayX,TXDelayZ);
toc

% tic;
% Recon5 = DASLUTV5(pL2,cRF3,pL2.RXDelay,TXDelayX,TXDelayZ);
% toc

tic; 
ReconC = beamform.computeDAScrfBMode(cRF); 
toc


tic;
RXDelayChkCPP = RXLUTTest(pL,pL2.RXDelay,pL.RXDelay);
toc


%% Plot

figure
plotGammaScaleImage(Recon,0.5);

figure
plotGammaScaleImage(Recon2,0.5);

figure
plotGammaScaleImage(Recon3,0.5);

figure
plotGammaScaleImage(Recon4,0.5);

% figure
% plotGammaScaleImage(Recon5,0.5);

figure
plotGammaScaleImage(ReconC,0.5);

%% LUT Testing
RXDelayV1 = reshape(pL.RXDelay,[pL.szZ,pL.szX,pL.numEl]);
RXDelayV2 = pL2.RXDelay;
RXDelayChk = zeros(pL.szZ,pL.szX,pL.numEl,'single');

tic;
for ix = 1:pL.szX
    for iz = 1:pL.szZ
        
        ixD = p.szX - ix;
        for ie = 1:pL.numEl
            iD = ie + ixD
            
%             iD = ie + p.szX - ix;
            
            s = RXDelayV2(iD,iz);
            RXDelayChk(iz,ix,ie) = s;
        end
    end
end
toc


%% Helper Functions
function [pO] = initParamsLUT(p)

    pO = p;
    % --- LUT settings (factorized TX + RX, stored in samples) ---
    % We store delays in *sample units* (0-based, for direct Eigen indexing).
    % Total sample index (0-based) in MEX will be:
    %   s = RXDelay(pix,el) + TXDelay(pix,ang) + tShift
    %
    % where tShift = -t0*fs (accounts for acquisition start time offset).

    % Enforce types (important for MEX side speed + predictable casting)
    pO.fs     = double(pO.fs);
    pO.c      = double(pO.c);
    pO.t0     = double(pO.t0);

    pO.numEl  = int32(pO.numEl);
    pO.na     = int32(pO.na);
    pO.szX    = int32(pO.szX);
    pO.szZ    = int32(pO.szZ);
    pO.szAcq  = int32(pO.szAcq);

    % Flatten pixel grid in MATLAB's column-major order:
    % pix = iz + (ix-1)*szZ  (1-based in MATLAB)
    % We'll keep it as an Npoints-by-1 list to match Eigen column-major mapping.
    x = single(pO.xCoord(:));          % [szX,1]
    z = single(pO.zCoord(:));          % [szZ,1]
    [X,Z] = meshgrid(x, z);            % both [szZ, szX]
    Xv = X(:);                          % [nPoints,1]
    Zv = Z(:);                          % [nPoints,1]

    % Convert common scale: fs/c
    fsOverC = single(pO.fs / pO.c);

    % --- TX plane-wave delay LUT in samples: [nPoints x na] ---
    % Plane wave propagation delay: (x*sin(theta) + z*cos(theta)) / c
    % Convert to samples: fs/c * (...)
    theta = single(pO.TXangle(:)).';    % [1,na]
    pO.TXDelay = fsOverC * ((sign(theta)*p.L/2 + Xv) .* sin(theta) + Zv .* cos(theta));  % [nPoints,na]
    pO.TXDelay = single(pO.TXDelay) + single(pO.t0*pO.fs);

    % --- RX delay LUT in samples: [nPoints x numEl] ---
    % Spherical receive path: sqrt((x-x_el)^2 + z^2) / c
    % Convert to samples: fs/c * sqrt(...)
    el = single(pO.ElemPos(:)).';       % [1,numEl]
    dx = Xv - el;                        % implicit expansion -> [nPoints,numEl]
    pO.RXDelay = fsOverC * sqrt(dx.^2 + (Zv.^2));  % [nPoints,numEl]
    pO.RXDelay = single(pO.RXDelay);

    % Acquisition time-zero shift in samples (0-based indexing on MEX side)
%     pO.tShift = single(pO.t0 * pO.fs);
    
    
end


function [pO] = initParamsLUTV2(p)

    pO = p;
    % --- LUT settings (factorized TX + RX, stored in samples) ---
    % We store delays in *sample units* (0-based, for direct Eigen indexing).
    % Total sample index (0-based) in MEX will be:
    %   s = RXDelay(pix,el) + TXDelay(pix,ang) + tShift
    %
    % where tShift = -t0*fs (accounts for acquisition start time offset).

    % Enforce types (important for MEX side speed + predictable casting)
    pO.fs     = double(pO.fs);
    pO.c      = double(pO.c);
    pO.t0     = double(pO.t0);

    pO.numEl  = int32(pO.numEl);
    pO.na     = int32(pO.na);
    pO.szX    = int32(pO.szX);
    pO.szZ    = int32(pO.szZ);
    pO.szAcq  = int32(pO.szAcq);

    % Flatten pixel grid in MATLAB's column-major order:
    % pix = iz + (ix-1)*szZ  (1-based in MATLAB)
    % We'll keep it as an Npoints-by-1 list to match Eigen column-major mapping.
    x = single(pO.xCoord(:));          % [szX,1]
    z = single(pO.zCoord(:));          % [szZ,1]
    

    % Convert common scale: fs/c
    fsOverC = single(pO.fs / pO.c);

    % --- TX plane-wave delay LUT in samples: [nPoints x na] ---
    % Plane wave propagation delay: (x*sin(theta) + z*cos(theta)) / c
    % Convert to samples: fs/c * (...)
    theta = single(pO.TXangle(:)).';    % [1,na]
    pO.TXDelayX = (((sign(theta)*p.L/2 + x).*sin(theta))*fsOverC).';
    pO.TXDelayZ = (z.*cos(theta)*fsOverC + p.t0*p.fs).';

    % --- RX delay LUT in samples: [nPoints x numEl] ---
    % Spherical receive path: sqrt((x-x_el)^2 + z^2) / c
    % Convert to samples: fs/c * sqrt(...)
    el = single(pO.ElemPos(:)).';       % [1,numEl]
    
    pO.RXDelayZ = (z.^2);
    pO.RXDelayX = ((x-el).^2).';
    RXDelayX = pO.RXDelayX;
    
    % RX delay LUT with reduced redundancy and precomputed sqrt
    [m,n] = size(RXDelayX);
    ks = -(m-1):(n-1);

    repVals = zeros(length(ks),1);
    
    for idx = 1:numel(ks)
        d = diag(RXDelayX, ks(idx));
        repVals(idx) = mean(d);
    end
    
    pO.RXDelay = single(sqrt(repVals(:) + pO.RXDelayZ.')*fsOverC);
    
end
