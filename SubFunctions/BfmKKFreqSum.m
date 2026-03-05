function [Recon,RawDataKK] = BfmKKFreqSum(RFData, p)
%% BfmKKFreqSum - KK beamforming wrapper (frequency-domain compression)
%
% Purpose:
%   Performs KK beamforming by descrambling, FFT, frequency-domain
%   compression via CompressKKFourier, IFFT, and LUT-based beamforming
%   via KKLUT.
%
% Inputs:
%   RFData - raw RF data [samples x channels]
%   p      - parameter struct from initKKFreqDomain (must contain szAcq,
%            numEl, na, nRX, midpt, ConnMap, shiftFac, RXDelayX,
%            RXDelayZ, TXDelayX, TXDelayZ)
%
% Outputs:
%   Recon     - complex beamformed image [szZ x szX]
%   RawDataKK - compressed KK data after IFFT [szAcq x (na*nRX)]
%
% Pipeline:
%   1. Descramble: reorder channels via ConnMap, reshape to [szAcq x (na*numEl)]
%   2. FFT along time axis
%   3. Compress via CompressKKFourier MEX (frequency-domain shift-and-sum)
%   4. IFFT to return to time domain
%   5. Beamform via KKLUT MEX (plane-wave delay-and-sum over TX/RX angles)
%
% Dependencies: CompressKKFourier (MEX), KKLUT (MEX)
%
% See also: initKKFreqDomain, CompressKKFourier, KKLUT, BfmDAS

    % Descramble: reorder channels and reshape
    param = int32([p.szAcq,p.numEl,p.na,p.nRX,p.midpt]);
    RFData2 = reshape(RFData(1:(p.szAcq*p.na),p.ConnMap),[p.szAcq,p.na*p.numEl]);

    % FFT along time axis
    DataFFT = single(RFData2);
    DataFFT = fft(DataFFT,[],1);
    DataFFT = reshape(DataFFT,[p.szAcq,p.na*p.numEl]);

    % Frequency-domain compression (shift factors encode RX angle steering)
    DataInt = CompressKKFourier(param,DataFFT,p.shiftFac);

    % IFFT to return to time domain
    RawDataKK = ifft(DataInt,[],1);

    % LUT-based KK beamforming (plane-wave TX and RX delays)
    Recon = KKLUT(p,RawDataKK,p.RXDelayX,p.RXDelayZ,p.TXDelayX,p.TXDelayZ);

end
