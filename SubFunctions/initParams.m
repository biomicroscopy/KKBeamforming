function [param,RFdata,varargout] = initParams(dataFile,filetype)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   TODO: Nikunj will insert documentation here

    varargout = {};
    if (filetype == 0) 
        % Saved after exiting VSX with all data and structs in single file
        % input: cell containing filepath + filename in single entry
        % output: Scrambled angles RFData
        load(dataFile{1},'Trans','P','Resource','Receive','TW','TX','RcvData','PData');
        TXangle = reshape([TX(:).Steer],2,[]);
        TXangle = TXangle(1,:);
        na = length(TXangle);
        RFdata = RcvData{1}(:,:,1);
        
        VSonicsInit;
    elseif (filetype == 13) % init for boas lab data with L12.
        % Same as filetype 2, just accounts for the different aperture
        
        load(dataFile{1});
        
        TXangle = reshape([TX(:).Steer],2,[]);
        TXangle = TXangle(1,:);
        na = length(TXangle);
        RFdata = RcvData{1}(:,:,8);

        VSonicsInit;
        
        % Account for smaller aperture on tx/rx
        Apod = zeros(1,Trans.numelements);
        ap = TX(1).aperture;
        len = length(TX(1).Apod);
        Apod(ap:(ap+len-1)) = 1;
        Apod = logical(Apod);
        ElemPos = ElemPos(Apod);
        ConnMap = 1:len;
        numEl = len;

        RFdata = circshift(RFdata,-ap+1,2);
    end
    
    param = struct('fs',fs,... % sampling frequency [Hz]
            'pitch', pitch,... % Element Spacing [m]
            'fc', fc,... % center frequency [Hz]
            'c', c,... % speed of sound [m/s]
            'fnumber', fnumber,... % angular aperature ratio [ul]
            't0',t0,... % start of transmit [s]
            'TXangle',TXangle,... % vector of transmit angles [rad]
            'ElemPos',ElemPos,... % element position [m]
            'xCoord',xCoord,... % x-coordinates of grid [m]
            'zCoord',zCoord,... % z-coordinates of grid [m]
            'numEl',numEl,... % Number of elements [ul]
            'szRF',szRF,... %
            'szRFframe',szRFframe,... % Number of time samples in dataset - 1 [ul]
            'szX',szX,... % length of x-coordinate [ul]
            'szZ',szZ,... % length of z-coordinate [ul]
            'na',na,... % Number of beams [ul]
            'nc',nc,... % Number of channels [ul]
            'ConnMap',ConnMap,... % Element connector mapping [el]
            'startSample',startSample,... 
            'endSample',endSample,...
            'tShift',tShift,...
            'initFlag',1);
        
    param.nPoints = param.szX*param.szZ; % Number of pixels [ul]
    param.L = param.ElemPos(end)-param.ElemPos(1); % Total length of array [m]
    param.tShift = 0.0; 
    param.txPL = zeros(length(param.TXangle));
end

