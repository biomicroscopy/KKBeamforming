%% KKTimingExperiment - Benchmark KK vs DAS beamforming performance
%
% This script measures execution time of each pipeline stage (descramble,
% FFT, compress, IFFT, beamform) for KK and DAS beamforming across
% multiple RX angle counts (M) and TX angle counts. Runs 1000 iterations
% for statistical significance.
%
% User paths to modify:
%   - dataFilePath: path to ultrasound dataset directory
%
% Required data: ContrastTarget .mat files with 7 and 15 TX angles
%
% Required functions: initParams, initKKFreqDomain, initParamsLUTV2,
%   CompressKKFourier (MEX), KKLUT (MEX), DASLUT (MEX)
%
% Outputs: Table of mean execution times per pipeline stage

%% Initialize file location
clearvars

% Extract Current Path
currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

dataFile{1} = dataFilePath + "KK Data\TallPhantom_2.3.26\ContrastTarget_7A24R.mat";
filetype = 0;
[p1,RFData1] = initParams(dataFile,filetype);
p1.szAcq = int32(p1.szRFframe+1);
RFData1 = RFData1(1:p1.szAcq*p1.na,:);
p1.szRF = int32(p1.szAcq*p1.na);

dataFile{1} = dataFilePath + "KK Data\TallPhantom_2.3.26\ContrastTarget_15A24R.mat";
filetype = 0;
[p2,RFData2] = initParams(dataFile,filetype);
p2.szAcq = int32(p2.szRFframe+1);
RFData2 = RFData2(1:p2.szAcq*p2.na,:);
p2.szRF = int32(p2.szAcq*p2.na);
%% Initialize

n = 1;
Mvec = [7,19,21,57];
% Define RX angle confocal

for i = 1:length(Mvec)
    pK(n) = defineConfocalRXAng(p1,Mvec(i)); n = n+1;
end

for i = 1:length(Mvec)
    pK(n) = defineConfocalRXAng(p2,Mvec(i)); n = n+1;
end

tic;
p1DAS = initParamsLUTV2(p1);
p2DAS = initParamsLUTV2(p2);
toc
%% Speed Test
NTest = 1000;
T = zeros(NTest,5,length(Mvec)*2+2);

Too = tic;
for i = 1:NTest
    n = 1;
    
    % Small angle count KK
    for j = 1:length(Mvec)
        To = tic; 
        [Recon,T1,lab] = BfmKKFreqSumTimed(RFData1, pK(n)); 
        T(i,5,n) = toc(To);
        T(i,1:4,n) = T1;
        n = n+1;
    end
    
    % Small angle count DAS
    To = tic; 
    [Recon,T2,lab2] = BfmDASTimed(RFData1, p1DAS);
    T(i,5,n) = toc(To);
    T(i,1:4,n) = T2;
    n = n+1;
    
    % Large angle count KK
    for j = 1:length(Mvec)
        To = tic; 
        [Recon,T1,~] = BfmKKFreqSumTimed(RFData2, pK(n-1)); 
        T(i,5,n) = toc(To);
        T(i,1:4,n) = T1;
        n = n+1;
    end

    % Large angle count DAS
    To = tic; 
    [Recon,T2,~] = BfmDASTimed(RFData2, p2DAS);
    T(i,5,n) = toc(To);
    T(i,1:4,n) = T2;
    n = n+1;
    
end
toc(Too)

Tmean = squeeze(mean(T,1));

colNames = {'KK 7/7','KK 7/19','KK 7/21','KK 7/57','DAS 7','KK 15/7','KK 15/19','KK 15/21','KK 15/57','DAS 15'};
lab{1} = 'DescrambleAnd/OrFFT';
lab{2} = 'CompressAnd/OrHilbert';
TabulatedTimes = array2table(round(Tmean*1000,1),'RowNames',lab,'VariableNames',colNames)


%% Helper Functions
function [pK] = defineConfocalRXAng(p,M)

dTheta = mean(diff(p.TXangle));
o = (-floor(M/2):floor(M/2));
p.RXangle = sign(o).*dTheta.*(2*abs(o)/M + mod(abs(o),double(floor(p.na/2))));
p.nRX = int32(length(p.RXangle));

tic; pK = initKKFreqDomain(p); toc

end

function [Recon,Timing,lab] = BfmKKFreqSumTimed(RFData, p)

    n = 1;

    % Descramble
    tic;
    param = int32([p.szAcq,p.numEl,p.na,p.nRX,p.midpt]);
    RFData2 = reshape(RFData(1:(p.szAcq*p.na),p.ConnMap),[p.szAcq,p.na*p.numEl]);
    
    % FFT
    DataFFT = single(RFData2);
    DataFFT = fft(DataFFT,[],1);
    DataFFT = reshape(DataFFT,[p.szAcq,p.na*p.numEl]);
    Timing(n) = toc; n = n+1;

    % Compress and Hilb
    tic;
    DataInt = CompressKKFourier(param,DataFFT,p.shiftFac);
    Timing(n) = toc; n = n+1;
    
    % IFFT
    tic;
    RawDataKK = ifft(DataInt,[],1); 
    Timing(n) = toc; n = n+1;
    
    % Beamform
    tic;    
    Recon = KKLUT(p,RawDataKK,p.RXDelayX,p.RXDelayZ,p.TXDelayX,p.TXDelayZ);
    Timing(n) = toc; n = n+1;

    lab = {'DescrambleAndFFT','Compress','IFFT','Beamform','Total'};
end

function [Recon,Timing,lab] = BfmDASTimed(RFData, p)

    n = 1;

    % FFT
    tic;
    RFData2 = single(reshape(RFData,[p.szAcq,p.na*p.nc]));
    DataFFT = fft(RFData2,[],1);
    Timing(n) = toc; n = n+1;
    
    % Hilbert Xform
    tic;
    DataFFT = DataFFT.*p.h;
    Timing(n) = toc; n = n+1;
    
    % IFFT
    tic;
    Data = ifft(DataFFT,[],1);
    Data = reshape(Data,[p.szRF,p.nc]);
    Timing(n) = toc; n = n+1;
    
    % Beamform
    tic;
    Recon = DASLUT(p,Data,p.RXDelay,p.TXDelayX,p.TXDelayZ);
    Timing(n) = toc; n = n+1;
    
    lab = {'DescrambleAndFFT','Compress','IFFT','Beamform','Total'};


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
    pO = convertParamsToSingle(pO);
    pO = convertParamsToInteger(pO);

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
    pO.TXDelayX = (((sign(theta)*p.L/2 + x).*sin(theta))*fsOverC);
    pO.TXDelayZ = (z.*cos(theta)*fsOverC + p.t0*p.fs);

    % --- RX delay LUT in samples: [nPoints x numEl] ---
    % Spherical receive path: sqrt((x-x_el)^2 + z^2) / c
    % Convert to samples: fs/c * sqrt(...)
    el = single(pO.ElemPos(:)).';       % [1,numEl]
    
    pO.RXDelayZ = (z.^2)*(fsOverC^2);
    pO.RXDelayX = (((x-el).^2)*(fsOverC^2)).';
    RXDelayX = pO.RXDelayX;
    
    % RX delay LUT with reduced redundancy and precomputed sqrt
    [m,n] = size(RXDelayX);
    ks = -(m-1):(n-1);

    repVals = zeros(length(ks),1);
    
    for idx = 1:numel(ks)
        d = diag(RXDelayX, ks(idx));
        repVals(idx) = mean(d);
    end
    
    pO.RXDelay = single(sqrt(repVals(:) + pO.RXDelayZ.'));
    
    
    % Hilbert prefactor
    h = zeros(p.szAcq,1,'single');
    if p.szAcq > 0 && 2*fix(p.szAcq/2) == p.szAcq
        % even and nonempty
        h([1 p.szAcq/2+1]) = 1;
        h(2:p.szAcq/2) = 2;
    elseif p.szAcq>0
        % odd and nonempty
        h(1) = 1;
        h(2:(p.szAcq+1)/2) = 2;
    end
    pO.h = h;
    
end