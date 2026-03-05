function [pO] = initKKFreqDomain(p)
%% initKKFreqDomain - Initialize KK frequency-domain compression parameters
%
% Purpose:
%   Computes the shift factors (complex exponential phase ramps) and
%   delay lookup tables needed for frequency-domain KK beamforming.
%   The shift factors encode the RX angle steering and are applied in
%   the Fourier domain by CompressKKFourier. The TX and RX delay LUTs
%   use a plane-wave propagation model with separable lateral (X) and
%   axial (Z) components.
%
% Inputs:
%   p - parameter struct from initParams with additional fields:
%         nRX     - number of RX angles (int32)
%         RXangle - vector of receive steering angles [rad]
%
% Outputs:
%   pO - copy of p with additional fields:
%          shiftFac - complex single [(numEl*nRX) x szAcq]
%                     Precomputed shift factors with Hilbert weighting
%          midpt    - one-sided FFT midpoint index
%          TXDelayX - TX lateral delay [szX x na], in samples
%          TXDelayZ - TX axial delay   [szZ x na], in samples
%          RXDelayX - RX lateral delay [szX x nRX], in samples
%          RXDelayZ - RX axial delay   [szZ x nRX], in samples
%
% Algorithm:
%   Shift factors:
%     For each frequency bin k and element-RXangle pair (x, n):
%       shiftFac(x,n,k) = exp(-j*2*pi*k*nShifts(x,n)/N) * h(k)
%     where nShifts encodes the per-element phase slope for each RX
%     angle, and h(k) is the one-sided Hilbert spectral weight.
%   Delay LUTs:
%     TX and RX both use plane-wave model (separable sin/cos terms).
%
% Dependencies: convertParamsToSingle, convertParamsToInteger
%
% See also: initParams, CompressKKFourier, KKLUT, BfmKKFreqSum

    pO = p;
    pO.nRX = int32(p.nRX);

    % Precompute shift factor parameters
    s = double(2*p.pitch*p.fs/p.c);    % normalized element spacing in samples
    tSize = double(p.szAcq);            % FFT length
    xSize = double(p.numEl);             % number of elements
    RXangle = double(p.RXangle(:)).';    % [1 x nRX]

    % Frequency bin indices (centered, then ifftshifted for FFT ordering)
    k = ifftshift((-floor(tSize/2):ceil(tSize/2)-1).');
    k = permute(k(:),[3,2,1]);          % [1 x 1 x szAcq] for broadcasting

    % Per-element phase slopes for each RX angle
    slopes = s*sin(RXangle)/2;
    nShifts = abs(slopes).*abs((1:xSize).' - xSize.^((1-sign(slopes))/2));

    % One-sided Hilbert spectral weighting vector
    h = zeros(p.szAcq,1);
    if p.szAcq > 0 && 2*fix(p.szAcq/2) == p.szAcq
        % even length
        h([1 p.szAcq/2+1]) = 1;
        h(2:p.szAcq/2) = 2;
        pO.midpt = p.szAcq/2+1;
    elseif p.szAcq>0
        % odd length
        h(1) = 1;
        h(2:(p.szAcq+1)/2) = 2;
        pO.midpt = (p.szAcq+1)/2;
    end

    % Build shift factors: phase ramp * Hilbert weight
    % Result: [numEl x nRX x szAcq], then reshaped to [(numEl*nRX) x szAcq]
    shiftFac = exp(-1i*2*pi*k.*nShifts/tSize).*permute(h,[3,2,1]);
    shiftFac = single(shiftFac);
    pO.shiftFac = reshape(shiftFac,[p.numEl*p.nRX,p.szAcq]);

    % Enforce types for MEX compatibility
    pO = convertParamsToInteger(pO);
    pO = convertParamsToSingle(pO);

    x = single(pO.xCoord(:));          % [szX x 1]
    z = single(pO.zCoord(:));          % [szZ x 1]

    % Common scale factor: convert meters to samples
    fsOverC = single(pO.fs / pO.c);

    % --- TX plane-wave delay LUT (separable X and Z components) ---
    theta = single(pO.TXangle(:)).';    % [1 x na]
    pO.TXDelayX = (((sign(theta)*p.L/2 + x).*sin(theta))*fsOverC);
    pO.TXDelayZ = (z.*cos(theta)*fsOverC + p.t0*p.fs);

    % --- RX plane-wave delay LUT (separable X and Z components) ---
    thetaR = single(pO.RXangle(:)).';   % [1 x nRX]
    pO.RXDelayX = (((sign(thetaR)*p.L/2 + x).*sin(thetaR))*fsOverC);
    pO.RXDelayZ = (z.*cos(thetaR)*fsOverC);

end
