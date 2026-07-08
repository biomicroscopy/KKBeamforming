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

ipath = ['-I' fullfile(currentDir{1:nL-2},'inc','eigen-3.4.0')]; % mfilename also includes the filename in the directory, so we need to remove that
%% Compile
mingwFlags = {'CXXFLAGS="$CXXFLAGS -march=native -std=c++14 -fno-math-errno -ffast-math -fopenmp -DNDEBUG -w -Wno-error"',...
            'LDFLAGS="$LDFLAGS -fopenmp"','CXXOPTIMFLAGS="-O3"'};

tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'DASLUT.cpp'); toc

%% Load Data
currentDir = matlab.desktop.editor.getActiveFilename;
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"KKBeamforming"),1)},"Datasets");

dataFile{1} = fullfile(dataFilePath, "ResolutionTargets_48.mat");
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

cRF3 = single(cRF);
%% Execute

tic;
Recon = DASLUT(pL2,cRF3,pL2.RXDelayX,pL2.RXDelayZ,pL2.TXDelayX,pL2.TXDelayZ,pL2.RXApX,pL2.RXApZ);
toc

%% Plot

figure
plotGammaScaleImage(Recon,0.5);
