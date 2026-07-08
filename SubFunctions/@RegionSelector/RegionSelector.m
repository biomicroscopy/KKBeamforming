classdef RegionSelector < handle
    %RegionSelector Class that handles selection,mapping and definition of
    %   metric regions

    properties
        
        % Grid Properties
        xCoord(1,:) double {mustBeReal, mustBeFinite, mustBeVector}
        zCoord(1,:) double {mustBeReal, mustBeFinite, mustBeVector}
        szX(1,1) int32 {mustBeInteger}
        szZ(1,1) int32 {mustBeInteger}

        % Region Properties
        nReg(2,2) double {mustBeReal, mustBeFinite}     % [x1,z1;x2,z2] coordss of opposite corners
        sReg(2,2) double {mustBeReal, mustBeFinite}     % [x0,rx;z0,rz] center and radii of ellipse
        iReg(2,2) double {mustBeReal, mustBeFinite}     % [x1,z1;x2,z2] coords of opposite corners

        nRegC(2,2) double {mustBeReal, mustBeFinite}
        sRegC(2,2) double {mustBeReal, mustBeFinite}
        iRegC(2,2) double {mustBeReal, mustBeFinite}

        % Other Properties
        name string
    end

    methods
        function obj = RegionSelector(param,nReg,sReg,iReg,name)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            % Set Grid coords
            obj.xCoord = param.xCoord;
            obj.zCoord = param.zCoord;
            obj.szX = param.szX;
            obj.szZ = param.szZ;
            
            % Set region indices
            obj.nReg = nReg;
            obj.sReg = sReg;
            obj.iReg = iReg;

            obj.name = name;

            obj.selectCoord();

        end

        % Change Grid Function - updates parameters based on a new input
        % grid. Throws an error if region isn't fully contained by new grid
        function changeGrid(obj,pNew)

            % Find new point indices based on absolute coordinates and new
            % grid. Then update all parameters
            dx = mean(diff(pNew.xCoord));
            minx = min(pNew.xCoord(:));
            dz = mean(diff(pNew.zCoord));
            minz = min(pNew.zCoord(:));

            nRegT(1) = round((obj.nRegC(1)-minx)/dx+1);
            nRegT(2) = round((obj.nRegC(2)-minx)/dx+1);
            nRegT(3) = round((obj.nRegC(3)-minz)/dz+1);
            nRegT(4) = round((obj.nRegC(4)-minz)/dz+1);

            sRegT(1) = round((obj.sRegC(1)-minx)/dx+1);
            sRegT(2) = round((obj.sRegC(2)-minz)/dz+1);
            sRegT(3) = obj.sRegC(3)/dx;
            sRegT(4) = obj.sRegC(4)/dz;

            iRegT(1) = round((obj.iRegC(1)-minx)/dx+1);
            iRegT(2) = round((obj.iRegC(2)-minx)/dx+1);
            iRegT(3) = round((obj.iRegC(3)-minz)/dz+1);
            iRegT(4) = round((obj.iRegC(4)-minz)/dz+1);

            % Check if regions are contained within new grid
            % Note: If the signal region is out of bounds, then it is 
            % necessary that the image region is out of bounds. Therefore, 
            % it is sufficient to check the image region only
            if (any(nRegT(:) <= 0) || any(sRegT(:) <= 0) || any(iRegT(:) <= 0))
                error("One of the resulting region dims was negative. Check lower bounds");
            elseif (nRegT(2) > pNew.szX || nRegT(4) > pNew.szZ)
                error("One of the noise region params falls outside array bounds");
            elseif (iRegT(2) > pNew.szX || iRegT(4) > pNew.szZ)
                error("One of the image region params falls outside array bounds");
            end 

            obj.xCoord = pNew.xCoord;
            obj.zCoord = pNew.zCoord;
            obj.szX = pNew.szX;
            obj.szZ = pNew.szZ;

            obj.nReg = reshape(nRegT,2,2);
            obj.sReg = reshape(sRegT,2,2);
            obj.iReg = reshape(iRegT,2,2);

        end


        % Extract map functions
        function [map] = getiMap(obj,rshpFlag)
            map = obj.mapRect(obj.iReg(3),obj.iReg(4),obj.iReg(1),obj.iReg(2));
            if rshpFlag
                map = reshape(map,[obj.szZ,obj.szX]);
            end
        end

        function [map] = getsMap(obj,rshpFlag)
            map = obj.mapCirc(obj.sReg(1),obj.sReg(2),obj.sReg(4),obj.sReg(3));
            if rshpFlag
                map = reshape(map,[obj.szZ,obj.szX]);
            end
        end

        function [map] = getnMap(obj,rshpFlag)
            map = obj.mapRect(obj.nReg(3),obj.nReg(4),obj.nReg(1),obj.nReg(2));
            if rshpFlag
                map = reshape(map,[obj.szZ,obj.szX]);
            end
        end

        % Metric Functions
        function [C, CNR, GCNR] = computeContMetrics(obj,images)

            nImg = length(images);
            signalMap = obj.getsMap(1);
            noiseMap = obj.getnMap(1);
            imgMap = obj.getiMap(1);
            C = zeros(nImg,1); CNR = C; GCNR = C;

            for n = 1:nImg
                if ~isreal(images(n).data) % Take magnitude value if complex
                    img = abs(images(n).data);
                else
                    img = images(n).data;
                end
                img = img./max(img(:));

                mu_i=mean(img(signalMap)); mu_o=mean(img(noiseMap));
                v_i=var(img(signalMap)); v_o=var(img(noiseMap));

                C(n)=10*log10(mu_i./mu_o);
                CNR(n)=10*log10(abs(mu_i-mu_o)/sqrt(v_i+v_o));

                % Select subset of image to generate histogram for. Compute w.r.t.
                % magnitudes
                imgReg = img(imgMap);
                x=linspace(min(imgReg(:)),max(imgReg(:)),100);

                [pdf_i]=hist(img(signalMap),x);
                [pdf_o]=hist(img(noiseMap),x);

                OVL=sum(min([pdf_i./sum(pdf_i); pdf_o./sum(pdf_o)]));
                GCNR(n)= 1 - OVL;

            end

        end


        % Plotting Functions
        function plotImageMarked(obj,zoomFlag)

            if zoomFlag
                nReg2 = obj.nReg;
                sReg2 = obj.sReg;
                obj.nReg = obj.nReg - obj.iReg(1,:);
                obj.sReg(:,1) = obj.sReg(:,1) - obj.iReg(1,:).';
            end


            theta = 0 : 0.01 : 2*pi;
            x = obj.sReg(3) * cos(theta) + obj.sReg(1);
            y = obj.sReg(4) * sin(theta) + obj.sReg(2);

            nRectW = obj.nReg(2)-obj.nReg(1);
            nRectH = obj.nReg(4)-obj.nReg(3);


            if ~zoomFlag
                iRectW = obj.iReg(2)-obj.iReg(1);
                iRectH = obj.iReg(4)-obj.iReg(3);
                rectangle('Position',[obj.iReg(1)+0.5,obj.iReg(3)+0.5,iRectW,iRectH],'EdgeColor','g')
            end
            rectangle('Position',[obj.nReg(1)+0.5,obj.nReg(3)+0.5,nRectW,nRectH],'EdgeColor','r')
            plot(x, y, 'LineWidth', 1,'Color','b');

            if zoomFlag
                obj.nReg = nReg2;
                obj.sReg = sReg2;
            end

        end
        
        function plotImageMarkedCoords(obj,scale,zoomFlag)

            theta = 0 : 0.01 : 2*pi;
            x = (obj.sRegC(3) * cos(theta) + obj.sRegC(1))*scale;
            y = (obj.sRegC(4) * sin(theta) + obj.sRegC(2))*scale;

            nRectW = (obj.nRegC(2)-obj.nRegC(1))*scale;
            nRectH = (obj.nRegC(4)-obj.nRegC(3))*scale;


            dx = mean(diff(obj.xCoord*scale));
            dz = mean(diff(obj.zCoord*scale));
            if ~zoomFlag
                iRectW = (obj.iRegC(2)-obj.iRegC(1))*scale;
                iRectH = (obj.iRegC(4)-obj.iRegC(3))*scale;
                rectangle('Position',[obj.iRegC(1)*scale+0.5*dx,obj.iRegC(3)*scale+0.5*dz,iRectW,iRectH],'EdgeColor','g')
            end
            rectangle('Position',[obj.nRegC(1)*scale+0.5*dx,obj.nRegC(3)*scale+0.5*dz,nRectW,nRectH],'EdgeColor','r')
            plot(x, y, 'LineWidth', 0.5,'Color','b');

        end
        
        function [imgout] = imgIdx(obj,img)
            imgout = img(obj(1).iReg(3):obj(1).iReg(4),obj(1).iReg(1):obj(1).iReg(2));
        end
        
        function [xout] = xIdx(obj)
            xout = obj.xCoord(obj(1).iReg(1):obj(1).iReg(2));
        end
        
        function [zout] = zIdx(obj)
            zout = obj.zCoord(obj(1).iReg(3):obj(1).iReg(4));
        end

    end

    methods(Access=protected,Sealed)

        % Function to define coordinates of regions based on indices
        function selectCoord(obj)

            dx = mean(diff(obj.xCoord));
            dz = mean(diff(obj.zCoord));

            obj.nRegC(1) = obj.xCoord(obj.nReg(1));
            obj.nRegC(2) = obj.xCoord(obj.nReg(2));
            obj.nRegC(3) = obj.zCoord(obj.nReg(3));
            obj.nRegC(4) = obj.zCoord(obj.nReg(4));

            obj.sRegC(1) = obj.xCoord(obj.sReg(1));
            obj.sRegC(2) = obj.zCoord(obj.sReg(2));
            obj.sRegC(3) = dx*(obj.sReg(3));
            obj.sRegC(4) = dz*obj.sReg(4);

            obj.iRegC(1) = obj.xCoord(obj.iReg(1));
            obj.iRegC(2) = obj.xCoord(obj.iReg(2));
            obj.iRegC(3) = obj.zCoord(obj.iReg(3));
            obj.iRegC(4) = obj.zCoord(obj.iReg(4));

        end


        % Generate Map functions
        function [map] = mapRect(obj,z1,z2,x1,x2)

            map = zeros(obj.szZ,obj.szX,"logical");
            rows = 1:obj.szZ; cols = 1:obj.szX;
            map( (rows > z1) & (rows <= z2), (cols > x1) & (cols <= x2) ) = 1;
        end

        function [map] = mapCirc(obj,x0,z0,rz,rx)

            % Create a logical image of an ellipse with specified diameter, center,
            % and image size. First create the image.
            [columnsInImage, rowsInImage] = meshgrid(1:obj.szX, 1:obj.szZ);
            % Next create the circle in the image.
            map = ((rowsInImage - z0)/rz).^2 + ((columnsInImage - x0)/rx).^2 <= 1;
        end





    end
end