%% VSonicsInit - Extract imaging parameters from Verasonics workspace
%
% Purpose:
%   Initialization script (not a function) that extracts ultrasound
%   system parameters from Verasonics data structures loaded into the
%   workspace. Called internally by initParams after loading .mat data.
%
% Expected workspace variables (from Verasonics .mat file):
%   Resource - system configuration (speed of sound, channels)
%   Trans    - transducer definition (frequency, spacing, elements)
%   Receive  - receive acquisition settings (sample rate, start/end)
%   TW       - transmit waveform (peak timing)
%   P        - imaging parameters (start/end depth)
%   PData    - pixel grid definition (optional, created if missing)
%   TXangle  - transmit angle vector (set by initParams before calling)
%   na       - number of TX angles (set by initParams before calling)
%
% Sets workspace variables:
%   xCoord, zCoord  - pixel grid coordinates [m]
%   fs              - sampling frequency [Hz]
%   pitch           - element spacing [m]
%   fc              - center frequency [Hz]
%   c               - speed of sound [m/s]
%   fnumber         - f-number (angular aperture ratio)
%   t0              - transmit start time offset [s]
%   ElemPos         - element positions [m]
%   numEl           - number of elements [int32]
%   szRF            - total samples per receive buffer frame
%   szRFframe       - samples per acquisition
%   szX, szZ        - pixel grid dimensions
%   nc              - number of receive channels
%   ConnMap         - element-to-channel mapping
%   startSample, endSample - per-angle sample bounds
%   tShift, initFlag       - initialization state variables
%
% See also: initParams

lambda = Resource.Parameters.speedOfSound/(Trans.frequency*1e6);

% Build pixel grid if not already defined in the data file
if (~exist('PData'))
    PData(1).PDelta = [Trans.spacing, 0, 0.5];
    PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
    PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
    PData(1).Size(3) = 1;      % single image page
    PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth];
end

% Convert pixel grid from wavelengths to meters
xCoord = (PData(1).Origin(1) + (0:PData(1).Size(2)-1)*PData(1).PDelta(1))*lambda;
zCoord = (PData(1).Origin(3) + (0:PData(1).Size(1)-1)*PData(1).PDelta(3))*lambda;

wvlToM = Resource.Parameters.speedOfSound/(Trans.frequency*1e6);

% System parameters
fs = Receive(1).decimSampleRate*1e6;
pitch = Trans.spacingMm*1e-3;
fc = Trans.frequency*1e6;
c = Resource.Parameters.speedOfSound;

% Compute f-number from transmit angular aperture
if (exist('TXangle'))
    NA = sin(range(TXangle)/2);
    fnumber = cot(2*asin(NA));
else
    fnumber = 1/0.6;
end

% Timing and geometry
t0 = (-Receive(1).startDepth + TW(1).peak + 2*Trans.lensCorrection)/(Trans.frequency*1e6);
ElemPos = Trans.ElementPos(:,1).'*wvlToM;
numEl = int32(Trans.numelements);
szRF=Resource.RcvBuffer.rowsPerFrame;
szRFframe=Receive(1).endSample-Receive(1).startSample;
szX=length(xCoord);
szZ=length(zCoord);
nc=Resource.Parameters.numRcvChannels;

% Element-to-channel connector mapping
if isfield(Trans,'Connector')
    ConnMap=Trans.Connector.';
else
    ConnMap = double(1:numEl);
end

% Per-angle sample bounds
startSample=[Receive(1:na).startSample];
endSample=[Receive(1:na).endSample];

tShift = 0.0;
initFlag=1;
