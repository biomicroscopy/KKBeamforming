function [p2] = computeNewGrid(p,xLims,zLims,varargin)
    p2 = p;
    if nargin>3 % For case where length of array is given
        szX = varargin{1};
        szZ = varargin{2};
        % For case where indices of endpoints are given
        if (isequal(floor(xLims),xLims) && isequal(floor(zLims),zLims))
            p2.xCoord = linspace(p.xCoord(xLims(1)),p.xCoord(xLims(2)),szX);
            p2.zCoord = linspace(p.zCoord(zLims(1)),p.zCoord(zLims(2)),szZ);
        else % For case where coordinates of endpoints are given
            p2.xCoord = linspace(xLims(1),xLims(2),szX);
            p2.zCoord = linspace(zLims(1),zLims(2),szZ);
        end
        p2.szX = szX;
        p2.szZ = szZ;
    elseif (nargin == 3) % For case where only indices of endpoints are given
        p2.xCoord = p.xCoord(xLims(1):xLims(2));
        p2.zCoord = p.zCoord(zLims(1):zLims(2));
        p2.szX = length(p2.xCoord);
        p2.szZ = length(p2.zCoord);
    else
        error("Wrong number of inputs")
    end
    p2.nPoints = p2.szX*p2.szZ;
end