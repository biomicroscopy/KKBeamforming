function [pS] = convertParamsToInteger(p)

pS = p;

pS.numEl = int32(p.numEl);
pS.szRF = int32(p.szRF);
if (isfield(p,'szRFframe'))
    pS.szRFframe = int32(p.szRFframe);
    pS.szAcq = int32(p.szRFframe+1);
else
    pS.szAcq = int32(p.szAcq);
end
pS.szX = int32(p.szX);
pS.szZ = int32(p.szZ);
if (isfield(p,'szY'))
    pS.szY = int32(p.szY);
end

pS.na = int32(p.na);
pS.nc = int32(p.nc);
pS.ConnMap = int32(p.ConnMap);
pS.startSample = int32(p.startSample);
pS.endSample = int32(p.endSample);
pS.nPoints = int32(p.nPoints);



end