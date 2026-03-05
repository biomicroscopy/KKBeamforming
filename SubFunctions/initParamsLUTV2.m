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