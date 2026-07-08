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
%     pO.nShifts = nShifts;

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
