function [param,RFdata,varargout] = initParams(dataFile,filetype)
%% initParams - Load ultrasound RF data and initialize parameter struct
%
% Purpose:
%   Loads raw RF data from a Verasonics .mat file and constructs the
%   parameter struct used by all downstream beamforming functions.
%   Supports two acquisition formats (filetypes 0 and 13).
%
% Inputs:
%   dataFile - cell array {filepath} pointing to a Verasonics .mat file
%   filetype - integer selecting the data format:
%                0  = Standard VSX export (full aperture)
%                13 = BOAS lab L12 probe (sub-aperture acquisition)
%
% Outputs:
%   param  - struct containing all beamforming parameters (see fields below)
%   RFdata - raw RF data matrix [samples x channels]
%
% Parameter struct fields:
%   fs          - sampling frequency [Hz]
%   pitch       - element spacing [m]
%   fc          - center frequency [Hz]
%   c           - speed of sound [m/s]
%   fnumber     - angular aperture ratio [unitless]
%   t0          - start-of-transmit time offset [s]
%   TXangle     - vector of transmit steering angles [rad]
%   ElemPos     - element positions along the array [m]
%   xCoord      - lateral pixel coordinates [m]
%   zCoord      - axial pixel coordinates [m]
%   numEl       - number of active elements [int32]
%   szRF        - total samples in the receive buffer per frame
%   szRFframe   - samples per single acquisition (endSample - startSample)
%   szX         - number of lateral pixels
%   szZ         - number of axial pixels
%   na          - number of transmit angles
%   nc          - number of receive channels
%   ConnMap     - element-to-channel connector mapping
%   startSample - per-angle start sample indices
%   endSample   - per-angle end sample indices
%   nPoints     - total number of image pixels (szX * szZ)
%   L           - total array aperture length [m]
%
% Dependencies: VSonicsInit (script, called internally)
%
% See also: initParamsLUTV2, initKKFreqDomain, VSonicsInit

    varargout = {};
    if (filetype == 0)
        % Standard VSX export: all Verasonics structs saved in one file
        load(dataFile{1},'Trans','P','Resource','Receive','TW','TX','RcvData','PData');
        TXangle = reshape([TX(:).Steer],2,[]);
        TXangle = TXangle(1,:);
        na = length(TXangle);
        RFdata = RcvData{1}(:,:,1);

        VSonicsInit;
    elseif (filetype == 13)
        % BOAS lab L12 probe: sub-aperture acquisition requiring
        % aperture remapping and channel reordering

        load(dataFile{1});

        TXangle = reshape([TX(:).Steer],2,[]);
        TXangle = TXangle(1,:);
        na = length(TXangle);
        RFdata = RcvData{1}(:,:,8);

        VSonicsInit;

        % Remap aperture: zero-pad full array, then select active elements
        Apod = zeros(1,Trans.numelements);
        ap = TX(1).aperture;
        len = length(TX(1).Apod);
        Apod(ap:(ap+len-1)) = 1;
        Apod = logical(Apod);
        ElemPos = ElemPos(Apod);
        ConnMap = 1:len;
        numEl = len;

        % Shift RF data to align with sub-aperture element indices
        RFdata = circshift(RFdata,-ap+1,2);
    end

    param = struct('fs',fs,... % sampling frequency [Hz]
            'pitch', pitch,... % element spacing [m]
            'fc', fc,... % center frequency [Hz]
            'c', c,... % speed of sound [m/s]
            'fnumber', fnumber,... % angular aperture ratio [unitless]
            't0',t0,... % start of transmit [s]
            'TXangle',TXangle,... % vector of transmit angles [rad]
            'ElemPos',ElemPos,... % element positions [m]
            'xCoord',xCoord,... % lateral pixel coordinates [m]
            'zCoord',zCoord,... % axial pixel coordinates [m]
            'numEl',numEl,... % number of elements [int32]
            'szRF',szRF,... % total samples per frame
            'szRFframe',szRFframe,... % samples per acquisition
            'szX',szX,... % number of lateral pixels
            'szZ',szZ,... % number of axial pixels
            'na',na,... % number of TX angles
            'nc',nc,... % number of receive channels
            'ConnMap',ConnMap,... % element-to-channel mapping
            'startSample',startSample,... % per-angle start sample indices
            'endSample',endSample,... % per-angle end sample indices
            'tShift',tShift,... % time shift [s]
            'initFlag',1);

    param.nPoints = param.szX*param.szZ; % total number of image pixels
    param.L = param.ElemPos(end)-param.ElemPos(1); % total array aperture [m]
    param.tShift = 0.0;
    param.txPL = zeros(length(param.TXangle));
end
