function [images] = bfmAndProcessFreq(p,RFData,M)
%% bfmAndProcessFreq - Run DAS and KK beamforming with multiple RX configurations
%
% Purpose:
%   Beamforms the same RF dataset using DAS and several KK configurations,
%   returning a struct array of images. KK variants differ in how RX
%   angles are chosen relative to TX angles. Also computes coherent and
%   incoherent combinations of selected KK images.
%
% Inputs:
%   p      - parameter struct from initParams (with szAcq, szRFframe, etc.)
%   RFData - raw RF data [samples x channels]
%   M      - number of RX angles (must be odd)
%
% Outputs:
%   images - struct array with fields:
%              .name    - string label for this image variant
%              .data    - beamformed image [szZ x szX]
%              .RXangle - RX angles used [1 x nRX]
%            Entries (in order):
%              1: DAS (reference baseline)
%              2: KK RXEqual (RX angles = TX angles)
%              3: KK j=0 (offset parameter j=0)
%              4: KK j=3 (offset parameter j=3)
%              5: KK j=6 (offset parameter j=6)
%              6: KK confocal (wrapped offset pattern)
%              7: KK Coherent (coherent sum of j=0,3,6)
%              8: KK Incoherent (incoherent sum of j=0,3,6)
%
% Dependencies: initParamsLUTV2, initKKFreqDomain, BfmDAS, BfmKKFreqSum
%
% See also: BfmDAS, BfmKKFreqSum, initParamsLUTV2, initKKFreqDomain

    if ((round(M/2) == M/2))
        error("M must be odd")
    end

    images = repmat(struct('name',[],'data',[],'RXangle',[]),6,1);
    n = 1;

    k0 = 2*pi*p.fc/p.c;

    % --- 1. DAS (reference baseline) ---
    images(n).RXangle = asin((double((-p.numEl/2):(p.numEl/2 -1))*2*pi/p.L)/k0);
    pD = initParamsLUTV2(p);
    images(n).data = BfmDAS(RFData,pD);
    images(n).name = 'DAS'; n = n+1;

    % --- 2. KK with RX angles equal to TX angles ---
    p.nRX = p.na;
    p.RXangle = p.TXangle;
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK RXEqual"; n = n+1;

    % --- 3-5. KK with offset parameter j (controls RX angle spread) ---
    jStart = n;
    dTheta = mean(diff(p.TXangle));
    o = (-floor(M/2):floor(M/2));       % symmetric index vector

    % j=0: minimum spread
    j = 0;
    p.RXangle = sign(o).*dTheta.*(2*abs(o)/M + j);
    p.nRX = length(p.RXangle);
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK j=0"; n = n+1;

    % j=3: moderate spread
    j = 3;
    p.RXangle = sign(o).*dTheta.*(2*abs(o)/M + j);
    p.nRX = length(p.RXangle);
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK j=3"; n = n+1;

    % j=6: wide spread
    jEnd = n;
    j = 6;
    p.RXangle = sign(o).*dTheta.*(2*abs(o)/M + j);
    p.nRX = length(p.RXangle);
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK j=6"; n = n+1;

    % --- 6. Confocal KK (wrapped offset pattern) ---
    dTheta = mean(diff(p.TXangle));
    M2 = M*3;
    o2 = (-floor(M2/2):floor(M2/2));
    p.RXangle = sign(o2).*dTheta.*(2*abs(o2)/M2 + mod(abs(o2),double(floor(p.na/2))));
    p.nRX = length(p.RXangle);
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK confocal"; n = n+1;

    % --- 7. Coherent combination of j=0, j=3, j=6 ---
    images(n).RXangle = reshape([images(jStart:jEnd).RXangle].',M,jEnd-jStart+1);
    images(n).data = double(abs(sum(reshape([images(jStart:jEnd).data],p.szZ,p.szX,jEnd-jStart+1),3)).^2);
    images(n).name = "KK Coherent"; n = n+1;

    % --- 8. Incoherent combination of j=0, j=3, j=6 ---
    images(n).RXangle = reshape([images(jStart:jEnd).RXangle].',M,jEnd-jStart+1);
    images(n).data = double(sum(abs(reshape([images(jStart:jEnd).data],p.szZ,p.szX,jEnd-jStart+1)).^2,3));
    images(n).name = "KK Incoherent"; n = n+1;

end
