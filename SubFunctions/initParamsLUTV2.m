function [pO] = initParamsLUTV2(p)
%% initParamsLUTV2 - Compute delay lookup tables for DAS beamforming
%
% Purpose:
%   Creates factored TX and RX delay lookup tables (LUTs) in sample
%   units for use by the DASLUT MEX function. The TX delay is split
%   into separable lateral (X) and axial (Z) components. The RX delay
%   exploits Toeplitz-like redundancy across elements and lateral
%   positions, reducing storage to a single 1-D table per axial depth.
%   Also computes the one-sided Hilbert transform spectral weighting
%   vector used to form the analytic signal prior to beamforming.
%
% Inputs:
%   p - parameter struct from initParams (must contain fs, c, t0,
%       TXangle, ElemPos, xCoord, zCoord, szAcq, L, numEl, etc.)
%
% Outputs:
%   pO - copy of p with additional fields:
%          TXDelayX - TX lateral delay [szX x na], in samples
%          TXDelayZ - TX axial delay   [szZ x na], in samples
%          RXDelay  - RX delay table   [(numEl+szX-1) x szZ], in samples
%          h        - Hilbert spectral weights [szAcq x 1]
%
% Algorithm:
%   TX delay (plane wave):
%     TXDelayX = (sign(theta)*L/2 + x) * sin(theta) * fs/c
%     TXDelayZ = z * cos(theta) * fs/c + t0*fs
%   RX delay (spherical, reduced redundancy):
%     Full table: RXDelayX(el,x) = ((x - x_el)^2) * (fs/c)^2
%     Exploit constant anti-diagonals (Toeplitz structure) to extract
%     a single representative value per offset, then combine with
%     RXDelayZ = z^2 * (fs/c)^2 and take sqrt.
%
% Dependencies: convertParamsToSingle, convertParamsToInteger
%
% See also: initParams, DASLUT, BfmDAS

    pO = p;

    % Enforce types for MEX compatibility
    pO = convertParamsToSingle(pO);
    pO = convertParamsToInteger(pO);

    x = single(pO.xCoord(:));          % [szX x 1]
    z = single(pO.zCoord(:));          % [szZ x 1]

    % Common scale factor: convert meters to samples
    fsOverC = single(pO.fs / pO.c);

    % --- TX plane-wave delay LUT (separable X and Z components) ---
    % Plane wave delay: (sign(theta)*L/2 + x)*sin(theta) * fs/c  [lateral]
    %                    z*cos(theta) * fs/c + t0*fs              [axial]
    theta = single(pO.TXangle(:)).';    % [1 x na]
    pO.TXDelayX = (((sign(theta)*p.L/2 + x).*sin(theta))*fsOverC);
    pO.TXDelayZ = (z.*cos(theta)*fsOverC + p.t0*p.fs);

    % --- RX delay LUT (reduced redundancy via Toeplitz structure) ---
    % Full squared-distance table: ((x - x_el)^2) * (fs/c)^2
    el = single(pO.ElemPos(:)).';       % [1 x numEl]
    pO.RXDelayZ = (z.^2)*(fsOverC^2);
    pO.RXDelayX = (((x-el).^2)*(fsOverC^2)).';
    RXDelayX = pO.RXDelayX;

    % Extract anti-diagonal means to collapse the Toeplitz-like table
    % into a vector of length (numEl + szX - 1)
    [m,n] = size(RXDelayX);
    ks = -(m-1):(n-1);

    repVals = zeros(length(ks),1);
    for idx = 1:numel(ks)
        d = diag(RXDelayX, ks(idx));
        repVals(idx) = mean(d);
    end

    % Combine lateral and axial components under the square root
    % Result: [(numEl+szX-1) x szZ] delay in samples
    pO.RXDelay = single(sqrt(repVals(:) + pO.RXDelayZ.'));

    % --- Hilbert transform spectral weighting vector ---
    % Multiplied with FFT to produce one-sided spectrum (analytic signal)
    h = zeros(p.szAcq,1,'single');
    if p.szAcq > 0 && 2*fix(p.szAcq/2) == p.szAcq
        % even length
        h([1 p.szAcq/2+1]) = 1;
        h(2:p.szAcq/2) = 2;
    elseif p.szAcq>0
        % odd length
        h(1) = 1;
        h(2:(p.szAcq+1)/2) = 2;
    end
    pO.h = h;

end
