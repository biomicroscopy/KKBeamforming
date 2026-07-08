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
    
    pO.RXDelayZ = (z.^2);
    pO.RXDelayX = ((x-el).^2).';
    RXDelayX = pO.RXDelayX;
    
    % RX delay LUT with reduced redundancy and precomputed sqrt. Only works
    % when x = el.
    [m,n] = size(RXDelayX);
    ks = -(m-1):(n-1);

    repVals = zeros(length(ks),1);
    
    for idx = 1:numel(ks)
        d = diag(RXDelayX, ks(idx));
        repVals(idx) = mean(d);
    end
    
    pO.RXDelay = single(sqrt(repVals(:) + pO.RXDelayZ.')*fsOverC);
    
    % --- f-number (dynamic receive aperture) gating, applied in the condensed LUT ---
    % row -> |x - el| = sqrt(repVals);  column -> depth z.  Keep where |x-el| <= z*0.5/fnumber.
    halfApSq = (z(:).' * (0.5 / pO.fnumber)).^2;   % [1 x szZ], squared half-aperture per depth
    pO.RXDelay(repVals(:) > halfApSq) = 0;         % [L x 1] > [1 x szZ] broadcast -> [L x szZ] mask
    
    pO.RXApX = abs(x(:).' - el(:));
    pO.RXApZ = z(:)*0.5/pO.fnumber;

    pO.fsOverC = single(fsOverC);

    
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
