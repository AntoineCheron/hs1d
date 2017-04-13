classdef SWATHobj
% Create swath profile object (SWATHobj)
%
% Syntax
%
%    SW = SWATHobj(DEM)
%    SW = SWATHobj(DEM,xy)
%    SW = SWATHobj(DEM,S)
%    SW = SWATHobj(DEM,'pn','pv'...)
%
%
% Description
%
%     SWATHobj creates a swath profile object, which can be used to
%     obtain statistics on terrain attributes such as elevation or slope
%     angles along a line or other directional x,y pairs.
%     
%     SWATHobj(DEM) opens a figure to interactively create a SWATHobj with 
%     the user defining a line with an arbitrary amount of nodes.
%
%     SWATHobj(DEM,xy) creates a SWATHobj from x,y data.
%
%     SWATHobj(DEM,S) creates a SWATHobj from a STREAMobj.
%
%     SWATHobj(DEM,'pn','pv'...) creates a SWATHobj and defines certain
%     parameter values that control the geometry of the SWATHobj.
%
%
% Input arguments
%
%     DEM    digital elevation model (Class: GRIDobj)
%     xy     directional x,y data pair (n x 2 array)
%     S      drainage network (Class: STREAMobj)
%     
% Parameter name/value pairs   {default}
%
%     'width'    scalar {1e4}
%            width of the swath profile in meters
%
%     'gap'    scalar {0}
%            width of a gap centered on the trace of the swath profile
%            within which no data is obtained
%
%     'dx'    scalar {cellsize of DEM}
%            resampling distance in the longitudinal direction of the swath 
%            profile. Provide empty matrix ([]) for no resampling at all.
%
%     'dy'    scalar {cellsize of DEM}
%            resampling distance in the transverse direction of the swath 
%            profile. Provide empty matrix ([]) for no resampling at all.
%
%     'keepnodes'    {false},true
%            switch to adjust resampling of points along swath profile to 
%            make sure the original nodes are included. If activated (true)
%            spacing of points along profile will most likely not be
%            unique.
%
%     'keepdist'    false,{true}
%            switch to adjust distance vector of swath profile to match
%            original distance vector of input data. If deactivated (false)
%            distances will change according to resampling and smoothing.
%
%     'keeptrace'    {false},true
%            switch to adjust trace of swath profile to match original
%            trace of input data. If activated (true) but trace of profile 
%            has been resampled and smoothed, the 
%
%     'smooth'    scalar {0}
%            optional smoothing of profile trace in x,y space. Number 
%            corresponds to the length of the filter in map units along the
%            swath profile. Numbers greater than dx will result in 
%            progressive smoothing up to the point that the filtfilt 
%            function, which is used for the smoothing, reports an error.
%
%     'smoothlongest'     false,{true}
%            If the filter length (parameter 'smooth') is too long with
%            respect to the swath profile's length, this parameter
%            determines if the smoothing is skipped (false) or performed
%            with the longest filter length possible (true).
%
%     'plot'    false,{true}
%            switch to plot some basic statistics along the swath profile,
%            calculated transverse to the swath profile. Statistics include
%            minimum, maximum value and arithmetic mean +/- standard
%            deviation.
%
% Output
%
%     SW     swath profile object (SWATHobj)
%
%
% SWATHobj properties:
%     
%     xy0       - raw x,y points used to create SWATHobj
%     zd0       - elevation and distance of raw points
%     dx        - longitudinal resampling interval (xy-unit of DEM)
%     dy        - transverse resampling interval (xy-unit of DEM)
%     width     - width of SWATHobj (xy-unit of DEM)
%     gap       - central gap along SWATHobj (xy-unit of DEM)
%     smooth    - cell with scalars corresponding to smoothing values
%     xy        - x,y points of central ine of SWATHobj
%     distx     - distance along SWATHobj
%     disty     - distance across SWATHobj
%     X         - X coorindates of SWATHobj data points
%     Y         - Y coorindates of SWATHobj data points
%     Z         - Z values of SWATHobj at data points
%     name      - name of SWATHobj
%     xyunit    - xyunit taken from DEM (GRIDobj)
%     zunit     - zunit taken from DEM (GRIDobj)
%     georef    - georeference structure taken from DEM (GRIDobj)
%
% 
% Example 1
%
%     % interactively create swath profile
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM)
%     
% Example 2
%
%     % create SWATHobj along a STREAMobj
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,A>100);
%     S = trunk(klargestconncomps(S,1));
%     SW = SWATHobj(DEM,S,'smooth',200);
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: 23. May, 2013

    
    properties
        xy0
        zd0
        dx
        dy
        width
        gap
        smooth
        xy
        distx
        disty
        X
        Y
        Z
        name
        xyunit
        zunit
        georef
    end
    
    methods
        
        function SW = SWATHobj(DEM,varargin)
            
            if mod(nargin,2)~=0 % SWATHobj is created interactively
                hfig = figure;
                imagesc(DEM), axis image
                [XY] = getline;
            else % SWATHobj is created from provided points or STREAMobj
                XY = varargin{1};
                varargin = varargin(2:end);
                hfig=[];
            end
            
            
            %% Parse inputs
            p = inputParser;
            p.FunctionName = 'SWATHobj';
            addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
            addParamValue(p,'width',1e3,@(x) isnumeric(x))
            addParamValue(p,'gap',0,@(x) isnumeric(x))
            addParamValue(p,'dx',DEM.cellsize,@(x) isnumeric(x))
            addParamValue(p,'dy',DEM.cellsize,@(x) isnumeric(x))
            addParamValue(p,'keepdist',true,@(x) islogical(x))
            addParamValue(p,'keeptrace',false,@(x) islogical(x))
            addParamValue(p,'keepnodes',false,@(x) islogical(x))
            addParamValue(p,'smooth',0,@(x) isnumeric(x))
            addParamValue(p,'smoothlongest',true,@(x) islogical(x))
            addParamValue(p,'plot',true,@(x) islogical(x))
            
            parse(p,DEM,varargin{:});

            % required
            DEM         = p.Results.DEM;
            % construct SWATHobj
            SW.width    = p.Results.width;
            SW.gap      = p.Results.gap;
            SW.dx       = p.Results.dx;
            SW.dy       = p.Results.dy;
            SW.georef   = DEM.georef;
            SW.name     = ['SWATHobj created from GRIDobj ',DEM.name];
            SW.zunit    = DEM.zunit;
            SW.xyunit   = DEM.xyunit;
            
            doplot      = p.Results.plot;
            keepdist    = p.Results.keepdist;
            keepnodes   = p.Results.keepnodes;
            keeptrace   = p.Results.keeptrace;
            smoothing   = p.Results.smooth;
            smoothlong  = p.Results.smoothlongest;
            
            
            %% Create line vector
            if isa(XY,'STREAMobj') % SWATHobj is created from STREAMobj
                validatealignment(XY,DEM);
                [x,y,d] = STREAMobj2XY(XY,XY.distance);
                x = flipud(x);
                y = flipud(y);
                d = flipud(d); 
                
            else % SWATHobj is created from x,y values
                x = XY(:,1);
                y = XY(:,2);
                d = getdistance(x,y);
            end
            
            
            %% Create SWATHobj
            nx = find(isnan(x));
            if isempty(nx)
                nx1 = 1;
                nx2 = length(x);
            else
                nx1 = nx+1;
                nx2 = [nx(2:end)-1;length(x)];
            end
            
            % Loop over line segements
            ct = 0;
            for i = 1 : length(nx1)
                x0 = x(nx1(i):nx2(i));
                y0 = y(nx1(i):nx2(i));
                d0 = d(nx1(i):nx2(i));
                
                if ~isempty(SW.dx) % resample along x
                    if keepnodes
                        this_x = []; this_y = []; this_d = [];
                        for k = 1 : length(x0)-1
                            [xt,yt,dt] = interpline(x0(k:k+1),y0(k:k+1),d0(k:k+1),SW.dx);
                            this_x = [this_x;xt];
                            this_y = [this_y;yt];
                            this_d = [this_d;d0(k)+dt];
                        end
                    else
                        [this_x,this_y,this_d] = interpline(x0,y0,d0,SW.dx);
                    end
                else % no resampling along x
                    this_x = x0;
                    this_y = y0;
                end
                ix0 = coord2ind(DEM,x0,y0);
                z0 = DEM.Z(ix0);
                
                if length(this_x)>1
                    ct = ct+1;
                    
                    SW.xy0{ct} = [x0,y0];
                    SW.zd0{ct} = [z0,d0];
                    if smoothing>0
                        sm = round(smoothing/SW.dx);
                        if sm>0
                            try
                                b = ones(1,sm)./sm;
                                this_xf = filtfilt(b,1,this_x);
                                this_yf = filtfilt(b,1,this_y);
                                smovalue = smoothing;
                            catch
                                if smoothlong
                                    fprintf(1,'Smoothing length reduced: not enough points.\n');
                                    sm = round((length(this_x)/3)-1);
                                    b = ones(1,sm)./sm;
                                    this_xf = filtfilt(b,1,this_x);
                                    this_yf = filtfilt(b,1,this_y); 
                                    smovalue = sm*SW.dx;
                                else
                                    fprintf(1,'Line smoothing skipped: not enough points.\n');
                                    this_xf = this_x;
                                    this_yf = this_y;
                                    smovalue = 0;
                                end
                            end
                            SW.smooth{ct} = smovalue;
                        else
                            fprintf(1,'Bad smoothing length: filter is zero.\n');
                        end
                    else
                        this_xf = this_x;
                        this_yf = this_y;
                    end
                    
                    dX = diff(this_xf); % dx between points along profile
                    dY = diff(this_yf); % dy between points along profile
                    
                    if ~keepdist
                        this_d = getdistance(this_xf,this_yf);
                    end
                    
                    if ~keeptrace
                        this_x = this_xf;
                        this_y = this_yf;
                    end
                    
                    SW.xy{ct} = [this_x,this_y];
                    SW.distx{ct} = this_d;
                    
                    hwidth = SW.width/2;
                    hgap = SW.gap/2;
                    if isempty(SW.dy)
                        SW.disty{ct} = -hwidth;
                    else
                        SW.disty{ct} = (-hwidth : SW.dy : -hgap)';
                    end
                    SW.disty{ct} = [SW.disty{ct}; flipud(abs(SW.disty{ct}))];
                    SW.disty{ct} = unique(SW.disty{ct},'stable');

                    % Create transverse profiles that are orthogonal to the center line
                    DX = [dX(1); dX];
                    DY = [dY(1); dY];
                    [theta,~] = cart2pol(DX,DY);  % direction between points along profile
                    theta_orthogonal = theta+pi/2; % orthogonals to center line of swath profile
                    [x_orthogonal,y_orthogonal] = pol2cart(theta_orthogonal,1); % dx, dy of orthogonals
                    % create new points
                    ny = length(SW.disty{ct});
                    nx = length(SW.distx{ct});
                    SW.X{ct} = repmat(SW.disty{ct},1,nx) .* repmat(x_orthogonal',ny,1) + repmat(this_x',ny,1);
                    SW.Y{ct} = repmat(SW.disty{ct},1,nx) .* repmat(y_orthogonal',ny,1) + repmat(this_y',ny,1);

                    % Interpolate GRIDobj values
                    SW.Z{ct} = nan(size(SW.X{ct}));
                    ix = 1:numel(SW.X{ct});
                    SW.Z{ct}(ix) = interp(DEM,SW.X{ct}(ix),SW.Y{ct}(ix));
                    
                end
            end
            
            if doplot
                if ishandle(hfig); hold on; 
                else figure; end
                % plot map view
                plot(SW), axis equal
                % plot profile view
                figure, plotdz(SW);
            end

            
        end
        
    end % methods
        
end % classdef
