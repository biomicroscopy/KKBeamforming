function [Recon] = BfmDAS(RFData, p)

    % FFT
    RFData2 = single(reshape(RFData,[p.szAcq,p.na*p.nc]));
    DataFFT = fft(RFData2,[],1);
    
    % Hilbert Xform
    DataFFT = DataFFT.*p.h;
    
    % IFFT
    Data = ifft(DataFFT,[],1);
    Data = reshape(Data,[p.szRF,p.nc]);
    
    % Beamform
    Recon = DASLUTV5(p,Data,p.RXDelay,p.TXDelayX,p.TXDelayZ);
end