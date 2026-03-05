function [Recon,RawDataKK] = BfmKKFreqSum(RFData, p)

    % Descramble
    param = int32([p.szAcq,p.numEl,p.na,p.nRX,p.midpt]);
    RFData2 = reshape(RFData(1:(p.szAcq*p.na),p.ConnMap),[p.szAcq,p.na*p.numEl]);
    
    % FFT
    DataFFT = single(RFData2);
    DataFFT = fft(DataFFT,[],1);
    DataFFT = reshape(DataFFT,[p.szAcq,p.na*p.numEl]);

    % MEX
    DataInt = CompressKKFourierV4(param,DataFFT,p.shiftFac);
    
    % IFFT
    RawDataKK = ifft(DataInt,[],1); 
    
    % Beamform
    Recon = KKLUTV3(p,RawDataKK,p.RXDelayX,p.RXDelayZ,p.TXDelayX,p.TXDelayZ);

end


