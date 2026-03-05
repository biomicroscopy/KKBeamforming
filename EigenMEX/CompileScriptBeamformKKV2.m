clearvars
close all


currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
nL = length(currentDir);
ipath = ['-I' fullfile(currentDir{1:nL-2},'Verasonics Integrated C Code','lib','eigen-3.4.0')]; % mfilename also includes the filename in the directory, so we need to remove that
% ipath2 = ['-I' fullfile(currentDir{1:nL-1},'inc')];
%% Compile
mingwFlags = {'CXXFLAGS="$CXXFLAGS -march=native -std=c++14 -fno-math-errno -ffast-math -fopenmp -DNDEBUG -w -Wno-error"',...
            'LDFLAGS="$LDFLAGS -fopenmp"','CXXOPTIMFLAGS="-O3"'};

% tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'BfKKCPPV3.cpp'); toc
% tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'BfKKCPPMTX.cpp'); toc
% 
tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'KKLUT.cpp'); toc
% tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'KKLUTV2.cpp'); toc
% tic; mex(ipath,mingwFlags{1},mingwFlags{2},mingwFlags{3},'KKLUTV3.cpp'); toc
%% Load
currentDir = matlab.desktop.editor.getActiveFilename; 
currentDir = regexp(currentDir, filesep, 'split');
dataFilePath = fullfile(currentDir{1:find(contains(currentDir,"Ultrasound"),1)},"Datasets\");

dataFile{1} = dataFilePath + "KK Data\TallPhantom_12.17.25\ResolutionTargets_48.mat";
filetype = 2;

[p,RFData] = initParams(dataFile,filetype);
p.szAcq = int32(p.szRFframe+1);
cRF = computeCRF(double(RFData),p);
%% Initialize
% M = double(p.na);
% dTheta = mean(diff(p.TXangle));
% o = (-floor(M/2):floor(M/2));
% % j = floor(p.na/2)-4;
% j = 0;
% p.RXangle = sign(o).*dTheta.*(2*abs(o)/M + j);
% p.nRX = int32(length(p.RXangle));

dTheta = mean(diff(p.TXangle));
M2 = 7*3;
o2 = (-floor(M2/2):floor(M2/2));
p.RXangle = sign(o2).*dTheta.*(2*abs(o2)/M2 + mod(abs(o2),double(floor(p.na/2))));
p.nRX = int32(length(p.RXangle));

tic; pK = initKKFreqDomain(p); toc

%% Beamform

tic; idxtMTX = BfmKKPixelCPP(cRF, p); toc
ReconKKO = reshape(squeeze(sum(idxtMTX,[1,2])),[p.szZ,p.szX]);

tic; idxtKK = BfmKKFreqMTX(RFData,pK); toc
ReconKK = reshape(squeeze(sum(idxtKK,[1,2])),[pK.szZ,pK.szX]);

T2 = tic; [ReconSum,~,TimingLUT,labLUT] = BfmKKFreqSumTimed(RFData,pK); TimingLUT(end+1) = toc(T2);
cell2struct(num2cell(TimingLUT),labLUT,2)

T = tic; [ReconSum,~,TimingLUT2,labLUT2] = BfmKKFreqSumTimedV2(RFData,pK); TimingLUT2(end+1) = toc(T);
cell2struct(num2cell(TimingLUT2),labLUT2,2)

T2 = tic; [ReconSum,~,TimingLUT,labLUT] = BfmKKFreqSumTimedV3(RFData,pK); TimingLUT(end+1) = toc(T2);
cell2struct(num2cell(TimingLUT),labLUT,2)





%% Plot

figure
plotGammaScaleImage(pK.xCoord*1e3,pK.zCoord*1e3,ReconKK,0.5);


figure
plotGammaScaleImage(pK.xCoord*1e3,pK.zCoord*1e3,ReconKKO,0.5);

figure
plotGammaScaleImage(pK.xCoord*1e3,pK.zCoord*1e3,ReconSum,0.5);


%% Timing
N = 1000;
T = zeros(N,6,2);

for i = 1:N
    
    T2 = tic; 
    [ReconSum,~,TimingLUT,labLUT] = BfmKKFreqSumTimed(RFData,pK); 
    T(i,6,1) = toc(T2);
    T(i,1:5,1) = TimingLUT;
    
    T2 = tic; 
    [ReconSum,~,TimingLUT,labLUT] = BfmKKFreqSumTimedV3(RFData,pK); 
    T(i,6,2) = toc(T2);
    T(i,1:5,2) = TimingLUT;
end

Tmean = squeeze(mean(T,1))

% Plot of beamforming time comparison. This shows relative stability of
% iterating over RFData first as opposed to iterating over pixels. Smaller
% chunks of data and delay tables can be passed to the inner loop leading
% to smaller demands on the cache, thus minimizing cahce misses. Average
% speed increases a bit, but more importantly, overall stability increases
% as well.
chk = squeeze(T(:,5,:));
figure
plot(chk)

%% Helper functions


function [Recon,RawDataKK,Timing,lab] = BfmKKFreqSumTimed(RFData, p)

    n = 1;

    % Descramble
    tic;
    param = int32([p.szAcq,p.numEl,p.na,p.nRX,p.midpt]);
    RFData2 = reshape(RFData(1:(p.szAcq*p.na),p.ConnMap),[p.szAcq,p.na*p.numEl]);
    Timing(n) = toc; n = n+1;
    
    % FFT
    tic;
    DataFFT = single(RFData2);
    DataFFT = fft(DataFFT,[],1);
    DataFFT = reshape(DataFFT,[p.szAcq,p.na*p.numEl]);
    Timing(n) = toc; n = n+1;

    % MEX
    tic;
    DataInt = CompressKKFourierV4(param,DataFFT,p.shiftFac);
    Timing(n) = toc; n = n+1;
    
    % IFFT
    tic;
    RawDataKK = ifft(DataInt,[],1); 
    Timing(n) = toc; n = n+1;
    
    % Beamform
    tic;
    Recon = KKLUT(p,RawDataKK,p.RXDelayX,p.RXDelayZ,p.TXDelayX,p.TXDelayZ);
    Timing(n) = toc; n = n+1;

    lab = {'DescrambleLUT','FFTLUT','CompressLUT','IFFTLUT','BeamformLUT','TotalLUT'};
end


function [Recon,RawDataKK,Timing,lab] = BfmKKFreqSumTimedV3(RFData, p)

    n = 1;

    % Descramble
    tic;
    param = int32([p.szAcq,p.numEl,p.na,p.nRX,p.midpt]);
    RFData2 = reshape(RFData(1:(p.szAcq*p.na),p.ConnMap),[p.szAcq,p.na*p.numEl]);
    Timing(n) = toc; n = n+1;
    
    % FFT
    tic;
    DataFFT = single(RFData2);
    DataFFT = fft(DataFFT,[],1);
    DataFFT = reshape(DataFFT,[p.szAcq,p.na*p.numEl]);
    Timing(n) = toc; n = n+1;

    % MEX
    tic;
    DataInt = CompressKKFourierV4(param,DataFFT,p.shiftFac);
    Timing(n) = toc; n = n+1;
    
    % IFFT
    tic;
    RawDataKK = ifft(DataInt,[],1); 
    Timing(n) = toc; n = n+1;
    
    TXDelayX = p.TXDelayX.';
    TXDelayZ = p.TXDelayZ.';
    RXDelayX = p.RXDelayX.';
    RXDelayZ = p.RXDelayZ.';
    
    % Beamform
    tic;
    Recon = KKLUTV3(p,RawDataKK,RXDelayX,RXDelayZ,TXDelayX,TXDelayZ);
    Timing(n) = toc; n = n+1;

    lab = {'DescrambleLUT','FFTLUT','CompressLUT','IFFTLUT','BeamformLUT','TotalLUT'};
end

function [Recon,RawDataKK,Timing,lab] = BfmKKFreqSumTimedV2(RFData, p)

    n = 1;

    % Descramble
    tic;
    param = int32([p.szAcq,p.numEl,p.na,p.nRX,p.midpt]);
    RFData2 = reshape(RFData(1:(p.szAcq*p.na),p.ConnMap),[p.szAcq,p.na*p.numEl]);
    Timing(n) = toc; n = n+1;
    
    % FFT
    tic;
    DataFFT = single(RFData2);
    DataFFT = fft(DataFFT,[],1);
    DataFFT = reshape(DataFFT,[p.szAcq,p.na*p.numEl]);
    Timing(n) = toc; n = n+1;
    
    % MEX
    tic;
    DataInt = CompressKKFourierV4(param,DataFFT,p.shiftFac);
    Timing(n) = toc; n = n+1;
    
    % IFFT
    tic;
    RawDataKK = ifft(DataInt,[],1); 
    Timing(n) = toc; n = n+1;
    
    % Beamform
    tic;
    RawDataKK = reshape(RawDataKK,[p.szAcq,p.na,p.nRX]);
    RawDataKK = permute(RawDataKK,[1,3,2]);
    RawDataKK = reshape(RawDataKK,[p.szAcq,p.nRX*p.na]);
    Timing(n) = toc; n = n+1;
    
    tic;
    Recon = KKLUTV2(p,RawDataKK,p.RXDelayX,p.RXDelayZ,p.TXDelayX,p.TXDelayZ);
    Timing(n) = toc; n = n+1;
    
    lab = {'DescrambleLUTV2','FFTLUTV2','CompressLUTV2','IFFTLUTV2','Permute','BeamformLUTV2','TotalLUTV2'};
end