function [h, udata] = scalebarY(varargin)

%Creates a scalebar on a set of axes
%EAS modified from Mathworks website 6.25.16
%commented out line 116: if this is a 2D plot, 3 vals are not needed for DataAspectRatio

% SCALEBAR
% SCALEBAR OFF
% SCALEBAR(PARAMETER, VALUE, ...)
% SCALEBAR(HAXES, PARAMETER, VALUE, ...)
% H = SCALEBAR(...)
% Draws a scalebar on the axes and returns handle to the scalebar. The
% DataAspectRatio property of the axes must be set to [1 1 1] and the view
% must be in 2D. All parameters are optional (note the default values
% below). SCALEBAR OFF deletes the current scalebar.
% 
% PARAMETER/VALUE pairs
%     hAxes:        handle to the axes (defaults to current axes)
%     ScaleLength:  length to show (in data units) (defaults to ~10% of the 
%                       x-axis limit range)
%     ScaleLengthRatio: ScaleLength/range(XLim)
%     Location:     location of the scalebar. Possible values are
%                       northeast (default)
%                       northwest
%                       southeast
%                       southwest
%                       [x y] data coordinates
%     Colour:       colour of scalebar in 1x3 RGB array (default is [0 0 0])
%     Bold:         draw with bold text and linewidth=2. 
%                       True or false(default)
%     Unit:         string containing units e.g. 'mm'
%
% Note: SCALEBAR sets the XLimMode and YLimMode of the axes to manual.

% Constants:
directions = {'northwest','northeast','southeast','southwest'};

% Set parameters to defaults
hAxes = gca;
%     scalelengthratio = 0.1;
%     scalelength = 100;                          % this will correspond to 1 ms (given that the orignal ss file had x values *10 = uSec)
%     direction = 'location';
    colour = [0 0 0];
    boldflag = false;
    linewidth = 2;
    fontweight = 'normal';
    unitstring = '';
%     DataAspectRatio = [1 1 1];

% PROCESS ARGUMENTS
    if nargin>0
        args = {};
        % Process if arguments given as "hAxes, ..."
        if ishandle(varargin{1}) %hAxes
            args{length(args)+1}='hAxes';
            args{length(args)+1}=varargin{1};
            varargin = varargin(2:end);
        end
        args=[args varargin];
        % Process argument pairs
        for n=1:2:length(args)
            parameter = args{n};
            if n == length(args)
                error('Parameter ''%s'' is not followed by a value in argument list', parameter);
            end
            value = args{n+1};
            switch lower(parameter)
                case 'haxes'
                    if ~ishandle(value)
                        error 'HAXES is not a valid Axes handle'
                    elseif ~strcmpi(get(value,'type'),'axes')
                        error 'HAXES is not a valid Axes handle'
                    elseif strcmpi(get(value,'tag'),'colorbar')
                        error 'HAXES is a handle to a colorbar'
                    else
                        hAxes = value;
                    end
                case 'scalelength'
                    if ~isnumeric(value)
                        error 'SCALELENGTH must be a numeric value'
                    end
                    scalelength = value;
                case 'scalelengthratio'
                    if ~isnumeric(value)
                        error 'SCALELENGTHRATIO must be a numeric value'
                    end
                    scalelengthratio = value;
                case 'location'
                    if ~(numel(value)==2 && isnumeric(value)) && ...
                       isempty(strmatch(lower(value),directions,'exact'))
                        error 'unrecognised value for LOCATION'
                    end
                    location = value;
                case 'colour'
                    if numel(value)~=3 || ~isnumeric(value)
                        error 'COLOUR must be a 1x3 representation of an RGB colour'
                    end
                    colour = value;
                case 'bold'
                    if ischar(value) && strcmpi(value,'true') || ...
                       (islogical(value) || isnumeric(value)) && value
                            boldflag = true;
                            linewidth=2;
                            fontweight='bold';
                    end
                case 'unit'
                    if ~ischar(value) || isempty(value)
                        error('''Unit'' must be followed by a string')
                    end
                    unitstring = [' ' strtrim(value)];
                otherwise
                    error(['unrecognised parameter: ' parameter]);
            end
        end
    end
    
    % CHECK IF VIEW IS IN 2D
    [az el] = view(hAxes);
    if el~=90
        error 'The Axes must be in 2D view'
    end
    
    %GET IMAGE AND AXES DATA
    axeslims = [get(hAxes,'xlim')' get(hAxes,'ylim')'];
    axesdir = [1 1];
    if strcmpi(get(hAxes,'XDir'),'reverse')
        axeslims(:,1) = flipud(axeslims(:,1));
        axesdir(1) = -1;
    end
    if strcmpi(get(hAxes,'YDir'),'reverse')
        axeslims(:,2) = flipud(axeslims(:,2));
        axesdir(2) = -1;
    end

    % CALCULATE SCALELENGTH
    if scalelength==0
        sl = range(axeslims(:,1))*scalelengthratio;
        slorder = 10^floor(log10(sl));
        scalelength = round(sl/slorder)*slorder;    
    else
        scalelengthratio = scalelength/range(axeslims(:,1));
    end
    
    %SET UP POSITIONING
    if ischar(location)
        switch location
            case 'northeast'
                anchor = [axeslims(2,1) - axesdir(1)*range(axeslims(:,1))*0.05, ...
                          axeslims(2,2) - axesdir(2)*range(axeslims(:,2))*0.05];
%                 direction = 'southwest';
                  direction = 'northeast';
            case 'northwest'
                anchor = [axeslims(1,1) + axesdir(1)*range(axeslims(:,1))*0.05, ...
                          axeslims(2,2) - axesdir(2)*range(axeslims(:,2))*0.05];
%                 direction = 'southeast';
                  direction = 'northwest';
            case 'southwest'
                anchor = [axeslims(1,1) + axesdir(1)*range(axeslims(:,1))*0.05, ...
                          axeslims(1,2) + axesdir(2)*range(axeslims(:,2))*0.05];
%                 direction = 'northeast';
                  direction = 'southwest';
            case 'southeast'
                anchor = [axeslims(2,1) - axesdir(1)*range(axeslims(:,1))*0.05, ...
                          axeslims(1,2) + axesdir(2)*range(axeslims(:,2))*0.05];
%                 direction = 'northwest';
                  direction = 'southeast';
        end    
    else
        anchor = location;
        if location
            dirToCentre = min(axeslims)+range(axeslims)/2 - location.*axesdir;
            direction = directions{ceil((-1*atan2(dirToCentre(2),dirToCentre(1))+pi)/(2*pi)*4)};
        end
    end

    linepos = [anchor; anchor];
    if ~isempty(strfind(direction,'northeast'))
%         linepos(2,1) = linepos(2,1)+axesdir(1)*scalelength; %creating the y-scale
        linepos(1,2) = linepos(2,1)-axesdir(1)*scalelength;   %for northeast
    else
        linepos(2,2) = linepos(1,1)-axesdir(1)*scalelength;   %for southwest
    end
    
    if ~isempty(strfind(direction,'north'))
        textalignment = {'bottom', 'center'};
    else
        textalignment = {'top', 'center'};
    end
            
    % DRAW SCALEBAR
    set(gca,'xlimmode','manual','ylimmode','manual');
    hg = hggroup('tag','scalebar');
    line(linepos(:,1), linepos(:,2), 'color', colour, 'linewidth', linewidth, 'parent', hg);

    if nargout>0
        h = hg;
    end
    
    % SETUP DELETE CALLBACK
    set(hg,'DeleteFcn',@deleteScaleBar)
    
    % SETUP LISTENER TO RESET SCALEBAR ON CHANGE OF AXES LIMITS
    hL(1) = addlistener(hAxes,'YLim','PostSet',@(src,event) resetScaleBar(src,event,hg));

    % SET USERDATA
    udata = {'ScaleLengthRatio',scalelengthratio;...
             'AnchorRatio',[(anchor(1)-min(axeslims(:,1)))/range(axeslims(:,1)) (anchor(2)-min(axeslims(:,2)))/range(axeslims(:,2))];...
             'Colour',colour;...
             'Listeners',hL;...
             'Bold',boldflag};
    set(hg,'UserData',udata);    
    
    % CALLBACK FUNCTIONS
    function deleteScaleBar(src,event)
        udata = get(src,'UserData');
        delete(udata{strcmpi(udata(:,1),'Listeners'),2});

    function resetScaleBar(src,event,SB)
        udata = get(SB,'UserData');
        hAxes = get(SB,'parent');
        
        delete(SB);        
        
        axeslims = [get(hAxes,'xlim')' get(hAxes,'ylim')'];
        
        scalelengthratio = udata{strcmpi(udata(:,1),'ScaleLengthRatio'),2};
        anchorratio = udata{strcmpi(udata(:,1),'AnchorRatio'),2};
        location = [anchorratio(1)*range(axeslims(:,1))+axeslims(1,1) anchorratio(2)*range(axeslims(:,2))+axeslims(1,2)];
        colour = udata{strcmpi(udata(:,1),'Colour'),2};
        boldflag = udata{strcmpi(udata(:,1),'Bold'),2};
        
        scalebar(hAxes,'ScaleLengthRatio',scalelengthratio,'Location',location,'Colour',colour, 'Bold', boldflag);
        
%         udataY = udata;
%  end

