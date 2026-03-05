function [images] = bfmAndProcessFreq(p,RFData,M)

    if ((round(M/2) == M/2))
        error("M must be odd")
    end

    images = repmat(struct('name',[],'data',[],'RXangle',[]),6,1);
    n = 1;

    k0 = 2*pi*p.fc/p.c;

    % Perform DAS
    images(n).RXangle = asin((double((-p.numEl/2):(p.numEl/2 -1))*2*pi/p.L)/k0);
    pD = initParamsLUTV2(p);
    images(n).data = BfmDAS(RFData,pD);
    images(n).name = 'DAS'; n = n+1;

    % Perform Equal RXangle KK
    p.nRX = p.na;
    p.RXangle = p.TXangle;
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK RXEqual"; n = n+1;

    % RXangle j = 0
    jStart = n;
    dTheta = mean(diff(p.TXangle));
    o = (-floor(M/2):floor(M/2));
    j = 0;
    p.RXangle = sign(o).*dTheta.*(2*abs(o)/M + j);
    p.nRX = length(p.RXangle);
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK j=0"; n = n+1;

    % RXangle j = 3
    j = 3;
    p.RXangle = sign(o).*dTheta.*(2*abs(o)/M + j);
    p.nRX = length(p.RXangle);
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK j=3"; n = n+1;

    % RXangle j = 7
    jEnd = n;
    j = 6;
    p.RXangle = sign(o).*dTheta.*(2*abs(o)/M + j);
    p.nRX = length(p.RXangle);
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK j=6"; n = n+1;

    % Confocal KK
    dTheta = mean(diff(p.TXangle));
    M2 = M*3;
    o2 = (-floor(M2/2):floor(M2/2));
    p.RXangle = sign(o2).*dTheta.*(2*abs(o2)/M2 + mod(abs(o2),double(floor(p.na/2))));
    p.nRX = length(p.RXangle);
    images(n).RXangle = p.RXangle;
    pK = initKKFreqDomain(p);
    images(n).data = BfmKKFreqSum(RFData, pK);
    images(n).data = double(images(n).data);
    images(n).name = "KK confocal"; n = n+1;
    
    % Coherent Combination
    images(n).RXangle = reshape([images(jStart:jEnd).RXangle].',M,jEnd-jStart+1);
    images(n).data = double(abs(sum(reshape([images(jStart:jEnd).data],p.szZ,p.szX,jEnd-jStart+1),3)).^2);
    images(n).name = "KK Coherent"; n = n+1;
    
    % Incoherent Combination
    images(n).RXangle = reshape([images(jStart:jEnd).RXangle].',M,jEnd-jStart+1);
    images(n).data = double(sum(abs(reshape([images(jStart:jEnd).data],p.szZ,p.szX,jEnd-jStart+1)).^2,3));
    images(n).name = "KK Incoherent"; n = n+1;
    

    


end