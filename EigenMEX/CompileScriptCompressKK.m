%% CompileScriptCompressKK - Compile and validate KK compression MEX
%
% This script compiles CompressKKFourier.cpp (MEX) and validates
% the frequency-domain compression against a pure MATLAB reference.
%
% User paths to modify:
%   - Eigen include path (ipath): points to Eigen 3.4.0 headers
%   - Data file path (dataFilePath): points to ultrasound dataset directory
%
% Required functions: initParams, initKKFreqDomain,
%   convertParamsToSingle, convertParamsToInteger, plotGammaScaleImage
%
% Dependencies: MinGW/GCC compiler with OpenMP support

clearvars
% close all


currentDir = matlab.desktop.editor.getActiveFilename;
currentDir = regexp(currentDir, filesep, 'split');
nL = length(currentDir);

ipath = ['-I' fullfile(currentDir{1:nL-3},'Verasonics Integrated C Code','lib','eigen-3.4.0')]; % mfilename also includes the filename in the directory, so we need to remove that
%% Compile
mingwFlags = {'CXXFLAGS="$CXXFLAGS -march=native -std=c++14 -fno-math-errno -ffast-math -fopenmp -DNDEBUG -w -Wno-error"',...
            'LDFLAGS="$LDFLAGS -fopenmp"','CXXOPTIMFLAGS="-O3"'};

tic; mex(ipath,mingwFlags{:},'CompressKKFourier.cpp'); toc

%% Load and initialize data

currentDir = matlab.desktop.editor.getActiveFilename;
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

dataFile{1} = dataFilePath + "KK Data\TallPhantom_12.17.25\ResolutionTargets_48.mat";
filetype = 0;

[p,RFData] = initParams(dataFile,filetype);
p.szAcq = int32(p.szRFframe+1);

%% Define RXangle
M = double(p.na);
dTheta = mean(diff(p.TXangle));
o = (-floor(M/2):floor(M/2));
j = floor(p.na/2)-1;
p.RXangle = sign(o).*dTheta.*(2*abs(o)/M + j);
p.nRX = length(p.RXangle);

s = 2*p.pitch*p.fs/p.c;

pK = initKKFreqDomain(p);

%% Test Run

tic; RawDataKKV1 = DataCompressV1(RFData, pK); toc

tic; RawDataKKV4 = DataCompressVMatlabPrecompute(RFData, pK); toc


%%
% Build full shiftFac from base phasor + h, then compare to pK.shiftFac

tSize = double(p.szAcq);
xSize = double(p.numEl);
RXSize = double(p.nRX);

base = exp(-1i*2*pi*abs(s*sin(double(p.RXangle))/2).*abs((1:xSize).' - xSize.^((1-sign(s*sin(double(p.RXangle))/2))/2))/tSize);
base = single(base);
full = complex(zeros(xSize,RXSize,tSize,'single'));

full(:,:,1) = 1;                     % k=0
for k = 2:tSize
    full(:,:,k) = full(:,:,k-1).*base;
end

% apply the SAME h used in initKKFreqDomain
h = zeros(tSize,1,'single');
if mod(tSize,2)==0
    h([1 tSize/2+1]) = 1;
    h(2:tSize/2) = 2;
else
    h(1) = 1;
    h(2:(tSize+1)/2) = 2;
end
full = full .* reshape(h,1,1,[]);

full_flat = reshape(full, [xSize*RXSize, tSize]); % (xSize*RXSize) x tSize
full_flat(:,1473) = conj(full_flat(:,1473));

err = max(abs(full_flat(:) - pK.shiftFac(:)));
disp(err)


%% Helper Functions
function [RawDataKK] = DataCompressV1(RFData, p)

    % Descramble
    param = int32([p.szAcq,p.numEl,p.na,p.nRX,p.midpt]);
    RFData2 = reshape(RFData(1:(p.szAcq*p.na),p.ConnMap),[p.szAcq,p.na*p.numEl]);

    % FFT
    DataFFT = single(RFData2);
    DataFFT = fft(DataFFT,[],1);
    DataFFT = reshape(DataFFT,[p.szAcq,p.na*p.numEl]);

    % MEX
    DataInt = CompressKKFourier(param,DataFFT,p.shiftFac);

    % IFFT
    RawDataKK = ifft(DataInt,[],1);

end


function [RawDataKK] = DataCompressVMatlabPrecompute(RFData, p)

    RFData2 = reshape(RFData(1:(p.szAcq*p.na),p.ConnMap),[p.szAcq,p.na*p.numEl]);

    DataFFT = single(RFData2);
    DataFFT = fft(DataFFT,[],1);
    DataFFT = reshape(DataFFT,[p.szAcq,p.na*p.numEl]);

    DataInt = complex(zeros(p.szAcq,p.na*p.nRX,'single'));  % pre-make complex

    for i = 1:p.midpt
        chk = reshape(DataFFT(i,:), p.na, p.numEl);

        chk2 = reshape(p.shiftFac(:,i), p.numEl, p.nRX);

        DataInt(i,:) = reshape(chk*chk2,[1,p.na*p.nRX]);

    end

    RawDataKK = ifft(DataInt,[],1);

end


function [pO] = initKKFreqDomain(p)

    pO = p;
    pO.nRX = int32(p.nRX);

    % precompute

    s = double(2*p.pitch*p.fs/p.c);
    tSize = double(p.szAcq);
    xSize = double(p.numEl);
    RXangle = double(p.RXangle(:)).';


    k = ifftshift((-floor(tSize/2):ceil(tSize/2)-1).');
    k = permute(k(:),[3,2,1]);

    slopes = s*sin(RXangle)/2;
    nShifts = abs(slopes).*abs((1:xSize).' - xSize.^((1-sign(slopes))/2));

    h = zeros(p.szAcq,1);
    if p.szAcq > 0 && 2*fix(p.szAcq/2) == p.szAcq
        % even and nonempty
        h([1 p.szAcq/2+1]) = 1;
        h(2:p.szAcq/2) = 2;
        pO.midpt = p.szAcq/2+1;
    elseif p.szAcq>0
        % odd and nonempty
        h(1) = 1;
        h(2:(p.szAcq+1)/2) = 2;
        pO.midpt = (p.szAcq+1)/2;
    end


    % Exponential
    shiftFac = exp(-1i*2*pi*k.*nShifts/tSize).*permute(h,[3,2,1]);
    shiftFac = single(shiftFac);
    pO.shiftFac = reshape(shiftFac,[p.numEl*p.nRX,p.szAcq]); % [numEl*nRX, szAcq]

    % Further checks
    pO = convertParamsToInteger(pO);

    pO = convertParamsToSingle(pO);

    % Flatten pixel grid in MATLAB's column-major order:
    % pix = iz + (ix-1)*szZ  (1-based in MATLAB)
    % We'll keep it as an Npoints-by-1 list to match Eigen column-major mapping.
    x = single(pO.xCoord(:));          % [szX,1]
    z = single(pO.zCoord(:));          % [szZ,1]


    % Convert common scale: fs/c
    fsOverC = single(pO.fs / pO.c);

    % --- TX plane-wave delay LUT in samples: [szX or szZ x na] ---
    % Plane wave propagation delay: (x*sin(theta) + z*cos(theta)) / c
    % Convert to samples: fs/c * (...)
    theta = single(pO.TXangle(:)).';    % [1,na]
    pO.TXDelayX = (((sign(theta)*p.L/2 + x).*sin(theta))*fsOverC);
    pO.TXDelayZ = (z.*cos(theta)*fsOverC + p.t0*p.fs);


    % --- RX plane-wave delay LUT in samples: [szX or szZ x nRX] ---
    % Plane wave propagation delay: (x*sin(theta) + z*cos(theta)) / c
    % Convert to samples: fs/c * (...)
    thetaR = single(pO.RXangle(:)).';    % [1,na]
    pO.RXDelayX = (((sign(thetaR)*p.L/2 + x).*sin(thetaR))*fsOverC);
    pO.RXDelayZ = (z.*cos(thetaR)*fsOverC);


end
