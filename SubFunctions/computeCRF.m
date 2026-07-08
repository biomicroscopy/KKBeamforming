function [cRF] = computeCRF(RFData,p)

    cRF = zeros(size(RFData),class(RFData));
    for i = 1:p.na
        cRF(p.startSample(i):p.endSample(i),:,:) = hilbert(RFData(p.startSample(i):p.endSample(i),:,:));
    end
    cRF([p.startSample(1:p.na)],:,:) = 0;
end