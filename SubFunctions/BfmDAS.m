function [Recon] = BfmDAS(RFData, p)
%% BfmDAS - Delay-And-Sum beamforming wrapper
%
% Purpose:
%   Performs DAS beamforming by applying a Hilbert transform (via FFT)
%   to produce an analytic signal, then calling the DASLUT MEX function
%   for LUT-based delay-and-sum reconstruction.
%
% Inputs:
%   RFData - raw RF data [samples x channels]
%   p      - parameter struct from initParamsLUTV2 (must contain szAcq,
%            na, nc, szRF, h, RXDelay, TXDelayX, TXDelayZ)
%
% Outputs:
%   Recon - complex beamformed image [szZ x szX]
%
% Pipeline:
%   1. Reshape RF data to [szAcq x (na*nc)]
%   2. FFT along time axis
%   3. Apply Hilbert spectral weights (one-sided spectrum)
%   4. IFFT to obtain analytic signal
%   5. Reshape to [szRF x nc] and beamform via DASLUT MEX
%
% Dependencies: DASLUT (MEX)
%
% See also: initParamsLUTV2, DASLUT, BfmKKFreqSum

    % FFT along time dimension
    RFData2 = single(reshape(RFData(1:(p.szAcq*p.na),:),[p.szAcq,p.na*p.nc]));
    DataFFT = fft(RFData2,[],1);
    
    % Hilbert Xform
    DataFFT = DataFFT.*p.h;
    
    % IFFT
    Data = ifft(DataFFT,[],1);
    Data = reshape(Data,[p.szAcq*p.na,p.nc]);
    p.szRF = int32(p.szAcq*p.na);

    % LUT-based DAS beamforming
    Recon = DASLUT(p,Data,p.RXDelayX,p.RXDelayZ,p.TXDelayX,p.TXDelayZ,p.RXApX,p.RXApZ);
end
