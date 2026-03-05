%% CompileAndTestLUTDAS - Compile and validate DAS LUT beamformer
%
% This script compiles DASLUT.cpp (MEX) and validates the output
% against a reference DAS beamformer implementation.
%
% User paths to modify:
%   - Eigen include path (ipath): points to Eigen 3.4.0 headers
%   - Data file path (dataFilePath): points to ultrasound dataset directory
%
% Required functions: initParams, computeCRF, convertParamsToSingle,
%   convertParamsToInteger, initParamsLUTV2, UnscrambleRF,
%   reconraw.DASBModeOffline, plotGammaScaleImage
%
% Dependencies: MinGW/GCC compiler with OpenMP support

clearvars
close all


currentDir = matlab.desktop.editor.getActiveFilename;
currentDir = regexp(currentDir, filesep, 'split');
nL = length(currentDir);

ipath = ['-I' fullfile(currentDir{1:nL-2},'Verasonics Integrated C Code','lib','eigen-3.4.0')]; % mfilename also includes the filename in the directory, so we need to remove that
%% Compile
mingwFlags = {'CXXFLAGS="$CXXFLAGS -march=native -std=c++14 -fno-math-errno -ffast-math -fopenmp -DNDEBUG -w -Wno-error"',...
            'LDFLAGS="$LDFLAGS -fopenmp"','CXXOPTIMFLAGS="-O3"'};

tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'DASLUT.cpp'); toc

%% Load Data
currentDir = matlab.desktop.editor.getActiveFilename;
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

dataFile{1} = dataFilePath + "KK Data\TallPhantom_12.17.25\ResolutionTargets_48.mat";
filetype = 0;

% Load and initialize data and parameter struct
[p,RFData] = initParams(dataFile,filetype);
p.szAcq = int32(p.szRFframe+1);

% Reorganize data into [samples x channels x angles] and hilbert transform
cRF = computeCRF(double(RFData),p);

%% Initialize
pL2 = initParamsLUTV2(p);

pL2 = convertParamsToSingle(pL2);
pL2 = convertParamsToInteger(pL2);

beamform = reconraw.DASBModeOffline(p);

cRF3 = single(cRF);
%% Execute

TXDelayX = pL2.TXDelayX.';
TXDelayZ = pL2.TXDelayZ.';

tic;
Recon = DASLUT(pL2,cRF3,pL2.RXDelay,TXDelayX,TXDelayZ);
toc

tic;
ReconC = beamform.computeDAScrfBMode(cRF);
toc

%% Plot

figure
plotGammaScaleImage(Recon,0.5);

figure
plotGammaScaleImage(ReconC,0.5);

%% Helper Functions

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
