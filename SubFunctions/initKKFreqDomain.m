function [pO] = initKKFreqDomain(p)

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



% Eigen::VectorXf dTX = (p.TXangle.array().sign()*p.L/2 + XPSF)*p.TXangle.array().sin() + ZPSF*p.TXangle.array().cos();
% Eigen::VectorXf dRX = (p.RXangle.array().sign()*p.L/2 + XPSF)*p.RXangle.array().sin() + ZPSF*p.RXangle.array().cos();
% 
% idxtKK = ((dRX.array() + dTX.transpose().array())*p.fs/p.c + p.t0*p.fs).round().matrix().cast<int>();
% 
% idxtKK = ( p.RXangle.array().sign()*(ZPSF*p.RXangle.array().tan() - XPSF) - p.L/2 <= 0 ).replicate(1,p.na).select(idxtKK,0);
% idxtKK = ( idxtKK.array() >= 0 ).select(idxtKK,0);
% idxtKK = ( idxtKK.array() <= p.szAcq-1 ).select(idxtKK,0);