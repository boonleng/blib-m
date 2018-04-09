% BLIB-M - Boonleng's library of functions
%
% 7/1/2017 - Corrected syntaxes for Matlab R2017b
%
% 6/2/2013 - Added rainbowmap().
%          - Added accentmap().
%
% 2/3/2013 - Added dmap().
%
% 12/6/2011 - Added figbrightness().
%
% 11/3/2011 - Added classmap().
%
% 8/28/2011 - Added drgmap().
%           - Added dgbmap().
%           - Added dbrmap().
%           - Added drgmapinv().
%           - Added dgbmapinv().
%           - Added dbrmapinv().
%           - Added singlebar().
%           - Added hsvdemo() - interactive HSV demonstration.
%
% 4/19/2011 - Code clean up.
%
% 4/15/2011 - Fixed a bug in spidergrid(). Now axis limits respect the
%             given input axes instead of gca.
%
% 8/16/2009 - Code clean up.
%           - Replaced some depreciating function calls.
%           - Added zmapn().
%           - Added zmapnx().
%           - Added rgmapn().
%
% 3/7/2009  - Added czmapx().
%
% 4/18/2008 - Added rcmap().
%           - Added bjetmapx().
%
% 4/13/2008 - Updated czmap().
%           - Added rgmap2().
%           - Added spidergrid.
%
% ------------------------------ Figures -------------------------------
% BMOVEFIG     - Move figure/axes without changing the size
%                boonlib('bmovefig',HANDLES,[X_OFFSET Y_OFFSET])
%
% BMOVEFIGANIM - Move figure/axes without changing the size
%                boonlib('bmovefig',HANDLES,[X_OFFSET Y_OFFSET])
%
% BZOOM        - Scale up/down a figure.
%                boonlib('bzoom',FIGURE_HANDLE,ZOOM_FACTOR)
%
% BSIZEWINDOW  - Resize the figure/axes and keep upper left position
%                boonlib('bsizewindow',FIGURE_HANDLE,[WIDTH HEIGHT])
%
% ----------------------------- Colormaps ------------------------------
%
% RBMAP        - Red to White to Blue Map.  Good for symmetric +/- range
%                (e.g. radial velocity)
%                cmap = boonlib('rbmap',NUMBER_OF_COLORS)
%
% GOMAP        - Green to Yellow to Orange, symmetric.
%                cmap = boonlib('gomap',NUMBER_OF_COLORS)
%
% OGMAP        - Orange to White to Green, symmetric.
%                cmap = boonlib('ogmap',NUMBER_OF_COLORS)
%
% RGMAP        - Red to Gray to Green Map.  Common for radial velocity
%                cmap = boonlib('rgmap',NUMBER_OF_COLORS)
%
% RGMAP2       - Red to Gray to Green Map 2.  Similar to one in NCDC Viewer.
%                Best for odd number of bins.
%                cmap = boonlib('rgmap2',NUMBER_OF_COLORS)
%
% RGMAPN       - Red to Gray to Green Map similar to NCDC Viewer.
%                Best for ven number of bins.
%                cmap = boonlib('rgmapn',NUMBER_OF_COLORS)
%
% RGMAPINV     - Using the above map, flipped around the center, Gray stays.
%                Good for highlighting middle shade
%                cmap = boonlib('rgmapinv',NUMBER_OF_COLORS)
%
% RAPMAP       - Colormap that prints well on both color and grayscale.
%                Black --> Blue --> Magenta --> Red --> Yellow --> White
%                Generated using ideas from Rappaport's paper.
%                cmap = boonlib('rapmap',NUMBER_OF_COLORS)
%
% ASYMAP       - Asymmetric Blue to White to Red Map.
%                (e.g. radial velocity with more emphasis on negative)
%                cmap = boonlib('asymap',[MIN MAX])
%                cmap = boonlib('asymap',[MIN MAX MINCLIP MAXCLIP])
%
% ZMAP         - A fairly standard reflectivity 16-shade colormap for radar images
%                cmap = boonlib('zmap') creates a colormap for [0 75] dBZ; or
%                cmap = boonlib('zmap',[MIN MAX]) for [MIN MAX] dBZ
%
% ZMAP2        - Reflectivity colormap used by Plymouth State Weather Center
%                cmap = boonlib('zmap2') creates a colormap for [0 75] dBZ; or
%                cmap = boonlib('zmap2',[MIN MAX]) for [MIN MAX] dBZ
%
% ZMAP3        - Reflectivity colormap used by WDSS-II
%                cmap = boonlib('zmap3') creates a colormap for [0 75] dBZ; or
%                cmap = boonlib('zmap3',[MIN MAX]) for [MIN MAX] dBZ
%
% ZMAPN        - Reflectivity colormap used by NOAA Weather and Climate Toolkit
%                cmap = boonlib('zmapn') creates a colormap for [-2.5 77.5] dBZ;
%
% CZMAP        - Reflectivity colormap with continuous shades of
%                blue, green, yellow, red, magenta
%
% CZMAPX       - Reflecivity colormap extended with continous shades of
%                dark magenta, charcoal, blue, green, yellow, red, magenta
%
% CARBMAP      - Carbone map. Suitable for refractivity.
%                cmap = boonlib('carbmap',NUMBER_OF_COLORS)
%
% CZMAP        - Similar to zmap with flexible number of colors.
%                cmap = boonlib('czmap',NUMBER_OF_COLORS)
%                Requirement: Minimum NUMBER_OF_COLORS = 16
%
% T7MAP        - Dark blue to Green at 0.7, Magenta at 0.7+ to yellow then white
%                cmap = boonlib('t7map',NUMBER_OF_COLORS)
%
% T6MAP        - Dark blue to Green at 0.6, Magenta at 0.6+ to yellow then white
%                cmap = boonlib('t6map',NUMBER_OF_COLORS)
%
% RCMAP        - Red, Orange, Yellow, White, Green, Dark Cyan, Cyan.
%                cmap = boonlib('rcmap',NUMBER_OF_COLORS)
%
% BJETMAP      - A cross between RGMAP and JET
%                cmap = boonlib('bjetmap',NUMBER_OF_COLORS)
%
% DMAP         - A colormap suitable for differential reflectivity
%                cmap = boonlib('dmap')
%
% RAINBOWMAP   - A colormap that resembles rainbow colors
%                cmap = boonlib('rainbowmap', NUMBER_OF_COLORS)
%
% BJETMAPX     - Extended BJETMAP
%                cmap = boonlib('bjetmapx',NUMBER_OF_COLORS)
%
% MAPINHEX     - Conver colormap from 0..1 representations to 0 to 255 int8
%                in HEX representation
%                boonlib('mapinhex',COLOR_MAP)
%
% ----------------------------- Drawings -------------------------------
%
% DRAW3DBOX    - Draw a 3D Cube at (x y z) = (0 0 0) and with a size of [1 1 1].
%                h = boonlib('draw3dbox',[X_MIN X_MAX Y_MIN Y_MAX Z_MIN Z_MAX])
%
% DRAWARROW    - Draw a solid arrow using fill() function.
%                Position parameters are in normalized unit of the gcf figure.
%                h = boonlib('drawarrow',[X Y LENGTH ANGLE_RAD])
%                h = boonlib('drawarrow',[X Y LENGTH ANGLE_RAD THICKNESS])
%
% CIRCLE       - Draw circles
%                h = boonlib('circle',[X Y R],NPTS) draws concentric circles at
%                [X Y] with radii of R.  Each circle is drawn with NPTS points.
%                Note: [X Y Z] must be in the size of Nx3 for drawing N circles.
%                boonlib('circle',[X Y Z]) draws circles with default NPTS=90.
%
% CSPACE       - Demonstrate RGB and YUV components of a colormap
%                boonlib('cspace', COLORMAP)
%
% GRAYME       - Set colormap to a Gray version (Preview of Monochrome printout)
%                boonlib('grayme')
%                boonlib('grayme',FIGURE_HANDLE)
%
% SPIDERGRID   - Draws spider-like grids at specified ring spacing and
%                number of lines
%                boonlib('spidergrid')
%                boonlib('spidergrid',AXES,RING_SPACING,NUMBER_OF_LINES)
%
% ------------------------------- Misc ---------------------------------
%
% CHOOSEFILE   - Interface to choose file using dir()
%                boonlib('choosefile',FILEPATH,FILE_CRITERIA)
%                boonlib('choosefile',FILEPATH_FUNCTION,FILE_CRITERIA)
%
%

function [varargout] = blib(fn,varargin)
if ~exist('fn','var')
	eval('help boonlib')
    if isempty(whos), feval('demo'); return; end
else
    if nargout==0; feval(fn,varargin{:});
    elseif nargout==1, x = feval(fn,varargin{:}); varargout(1) = {x};
	elseif nargout==2, [x, y] = feval(fn,varargin{:}); varargout(1) = {x}; varargout(2) = {y};
    end
end
return

% ------------------------------- Demo ---------------------------------
function [x] = demo()
feval('mydefault')
h = feval('assignfig','demo');
figure(h)
cm = {'bjetmap', 'bjetmapinv', 'bjetmapxinv',...
    'carbonemap', 'zmap','zmap2', 'zmap3', ...
	'czmap', 'czmapx', 'czmapxinv', 'zmapn', 'zmapnx', 'bzmap',...
    'gomap', 'ogmap', 'refmap',...
    'rgmapinv', 'rgmap2', 'rgmapw', 'rgmapn',...
    't7map', 't6map', 't5map',...
    'bgraymap', 'wmap', 'tempmapinv',...
    'brmap', 'rcmap',...
    'asymap0', 'asymap',...
    'vrmap',...
    'drgmap', 'dgbmap', 'dbrmap', ...
    'drgmapinv', 'dgbmapinv', 'dbrmapinv', ...
	'classmap', ...
    'dmap', 'rainbowmap', ...
    'singlebar', 'hsvdemo', ...
	'figbrightness', ...
    };
idx = 1;
while (idx<length(cm))
    fprintf('Map: %s\n',cm{idx});
    cmap = feval(cm{idx});
    cspace(cmap);
    if (idx==1)
        bzoom(gcf,1.0); drawnow;
        bxwin();
    end
    drawnow;
    idx = idx+1;
end
clf
set(gcf,'Visible','Off')
bsizewin(gcf,[200 200])
cm = {'mapinhexq', 'mapinhex',...
    'draw3dbox','drawarrow','circle','spidergrid',...
    'grayme','greenme','sepiame'};
idx = 1;
while (idx<length(cm))
    feval(cm{idx});
    idx = idx+1;
end
close(gcf);
feval('aintersect',1);
feval('choosefile','.','randomefilenamethatwontexist');
x = feval('isoverlap',[0 0 200 200],[100 100 200 200]);
feval('mydefault2')
feval('mydefault3')
return;
    
% ------------------------------ Figures -------------------------------

function bmovefig(arg1,arg2)
if nargin<2, error('Must have two input.'); end
if size(arg2,1) ~= length(arg1), error('Dimension of arg1 and arg2 does not match.'); end
nfig = length(arg1);
tmp = zeros([nfig 4]); for iarg = 1:length(arg1), tmp(iarg,1:4) = get(arg1(iarg),'Position'); end
for iarg = 1:length(arg1)
    set(arg1(iarg),'Position',[tmp(iarg,1:2)+arg2(iarg,1:2) tmp(iarg,3:4)]);
end


function bmovefiganim(arg1,arg2)
if nargin<2, error('Must have two input.'); end
if size(arg2,1) ~= length(arg1), error('Dimension of arg1 and arg2 does not match.'); end
del_step = 5.95943165613589.\arg2;
nfig = length(arg1);
tmp = zeros([nfig 4]); for iarg = 1:length(arg1), tmp(iarg,1:4) = get(arg1(iarg),'Position'); end
for idx = cumsum(1-sin(0:0.5*pi/15:0.49*pi))
    for iarg = 1:length(arg1)
        set(arg1(iarg),'Position',[round(tmp(iarg,1:2)+idx*del_step(iarg,:)) tmp(iarg,3:4)]);
        pause(0.01);
    end
end
return


function bsizewin(h,newsize)
if nargin < 1; fprintf('Please put something.\n'); end
win_position = feval('get',h,'Position');
new_lower_right = win_position(1:2) + [0,win_position(4)-newsize(2)];
feval('set',h,'Position',[new_lower_right,newsize]);
return


function bzoom(h,scale)
if nargin==0; h = feval('gcf'); scale = 2; end
win_position = feval('get',h,'Position');
newsize = floor(win_position(3:4)*scale);
new_lower_right = win_position(1:2) + [0,win_position(4)-newsize(2)];
feval('set',h,'Position',[new_lower_right,newsize]);
return


function bxwin(h,cn,reorder)
if (nargin<1)||isempty(h); h = get(0,'Children'); end
scnsize = get(0,'ScreenSize');
% If there are two screen (==2690 pixels), use my office setting, otherwise single screen settings
if scnsize(3)>2560, x_off = 1280; if (nargin<2)||isempty(cn), cn = [50 80]; end
else
    x_off = 0;
	if (nargin<2)||isempty(cn)
        if (ispc)
            cn = [10 30];
        else
            cn = [10 30];
        end
	end
end
% Buffer space among windows
if (ispc)
    bf = [20 20];
else
    bf = [5 5];
end
% Y offset for the dock bar
y_off = 50;

% Only get figures that still exist
h = h(ishandle(h)&(h~=0));
nfig = numel(h);
if ~exist('reorder','var'), reorder = 1; end
if (reorder)
	% Sort the figure so that they are from big to small
    hsiz = zeros(1,nfig);
	for ifig = 1:nfig
		tmp = get(h(ifig),'Position');
		hsiz(ifig) = tmp(3)*tmp(4);
	end
	[~, I] = sort(hsiz,2,'descend');
	h = h(I);
end
% 1-st available space is always at the upper-right corner
avspace = scnsize(3:4)-cn;
% Define it first
target = avspace;
f = zeros([nfig 4]);
for ifig = 1:nfig
    f(ifig,1:4) = get(h(ifig),'Position');
    % Number of available space
	nspace = size(avspace,1);
    % Fudge factor to measure remaining area
    ff = zeros([nspace 1]);
    % If there is menubar, penalize more
    if ~strcmpi(get(h(ifig),'Menubar'),'figure'), ybar = 21;
    else, ybar = 72;
    end
    % Go through the available space to see if current window fits without overlaps
    for ispace = 1:nspace
    	ff(ispace) = (avspace(ispace,1)-f(ifig,3)-x_off)*(avspace(ispace,2)-f(ifig,4)-y_off-ybar);
        tmp = (avspace(ispace,1)-f(ifig,3)     >=(target(1:ifig-1,1)-f(1:ifig-1,3)) & ...
               avspace(ispace,1)-f(ifig,3)     <=(target(1:ifig-1,1)              ) & ...
               avspace(ispace,2)-f(ifig,4)-ybar>=(target(1:ifig-1,2)-f(1:ifig-1,4)) & ...
               avspace(ispace,2)-f(ifig,4)-ybar<=(target(1:ifig-1,2)              )) | ...
              (avspace(ispace,1)-f(ifig,3)     >=(target(1:ifig-1,1)-f(1:ifig-1,3)) & ...
               avspace(ispace,1)-f(ifig,3)     <=(target(1:ifig-1,1)              ) & ...
               avspace(ispace,2)               >=(target(1:ifig-1,2)-f(1:ifig-1,4)) & ...
               avspace(ispace,2)               <=(target(1:ifig-1,2)              )) | ...
              (avspace(ispace,1)               >=(target(1:ifig-1,1)-f(1:ifig-1,3)) & ...
               avspace(ispace,1)               <=(target(1:ifig-1,1)              ) & ...
               avspace(ispace,2)-f(ifig,4)-ybar>=(target(1:ifig-1,2)-f(1:ifig-1,4)) & ...
               avspace(ispace,2)-f(ifig,4)-ybar<=(target(1:ifig-1,2)              )) | ...
              (avspace(ispace,1)               >=(target(1:ifig-1,1)-f(1:ifig-1,3)) & ...
               avspace(ispace,1)               <=(target(1:ifig-1,1)              ) & ...
               avspace(ispace,2)               >=(target(1:ifig-1,2)-f(1:ifig-1,4)) & ...
               avspace(ispace,2)               <=(target(1:ifig-1,2)              )) | ...
              (target(1:ifig-1,1)-f(1:ifig-1,3)>=(avspace(ispace,1)-f(ifig,3)     ) & ...
               target(1:ifig-1,1)-f(1:ifig-1,3)<=(avspace(ispace,1)               ) & ...
               target(1:ifig-1,2)-f(1:ifig-1,4)>=(avspace(ispace,2)-f(ifig,4)-ybar) & ...
               target(1:ifig-1,2)-f(1:ifig-1,4)<=(avspace(ispace,2)               )) | ...
              (target(1:ifig-1,1)-f(1:ifig-1,3)>=(avspace(ispace,1)-f(ifig,3)     ) & ...
               target(1:ifig-1,1)-f(1:ifig-1,3)<=(avspace(ispace,1)               ) & ...
               target(1:ifig-1,2)              >=(avspace(ispace,2)-f(ifig,4)-ybar) & ...
               target(1:ifig-1,2)              <=(avspace(ispace,2)               )) | ...
              (target(1:ifig-1,1)              >=(avspace(ispace,1)-f(ifig,3)     ) & ...
               target(1:ifig-1,1)              <=(avspace(ispace,1)               ) & ...
               target(1:ifig-1,2)-f(1:ifig-1,4)>=(avspace(ispace,2)-f(ifig,4)-ybar) & ...
               target(1:ifig-1,2)-f(1:ifig-1,4)<=(avspace(ispace,2)               )) | ...
              (target(1:ifig-1,1)              >=(avspace(ispace,1)-f(ifig,3)     ) & ...
               target(1:ifig-1,1)              <=(avspace(ispace,1)               ) & ...
               target(1:ifig-1,2)              >=(avspace(ispace,2)-f(ifig,4)-ybar) & ...
               target(1:ifig-1,2)              <=(avspace(ispace,2)               ));
        %if ~isempty(find(tmp,1)), ff(ispace) = -999999; end
        if any(tmp), ff(ispace) = -99999; end
    end
    % Check the position scores (debug only)
	%disp([avspace repmat([f(ifig,3) f(ifig,4)],[size(avspace,1) 1]) ff])

    % Sort the score based on x, y, ff
    savsp = [avspace ff];
    % Find out which space can fit it better
    [savsp, jj] = sortrows(savsp,[1 2 3]);
    %disp(savsp);
    ispace = jj(find(savsp(:,3)>0,1,'last'));
    if isempty(ispace), ispace = jj(end); end
    
    target(ifig,1:2) = avspace(ispace,:)-[0 ybar];
    %scnspace(max(target(ifig,2)-f(ifig,4)+1,1):target(ifig,2),...
    %         max(target(ifig,1)-f(ifig,3)+1,1):target(ifig,1)) = ifig;
    if ispace==nspace
        avspace = [avspace(1:ispace-1,:); ...
                   avspace(ispace,:)-[f(ifig,3)+bf(1) 0]; ...
                   avspace(ispace,:)-[0 f(ifig,4)+bf(2)+ybar]];
    else
        avspace = [avspace(1:ispace-1,:); ...
                   avspace(ispace,:)-[f(ifig,3)+bf(1) 0]; ...
                   avspace(ispace,:)-[0 f(ifig,4)+bf(2)+ybar]; ...
                   avspace(ispace+1:nspace,:)];
    end
end
offset = target-f(:,1:2)-f(:,3:4);
% figure(1)
% imagesc(scnspace);
% set(gca,'YDir','Normal')
if any(offset(:)~=0) 
    if verLessThan('matlab', '6.0.0') || (ispc)
        feval('bmovefig',h,offset);
    else
        feval('bmovefiganim',h,offset);
    end
end
return

% ----------------------------- Colormaps ------------------------------

function [x] = mybgcolor()
% tmp = get(0,'defaultTextColor');
% tmp = [0.3 0.59 0.11]*tmp(:);
%if tmp>0.5
if feval('figbrightness') < 0.5
	%x = [0.27 0.23 0.20];
	%x = [0.13 0.11 0.10];
	x = [0.0 0.0 0.0];
else
	% x = [0.925 0.900 0.850];
	x = [0.89 0.87 0.83];
	%x = [0.45 0.43 0.41];
end
return


function [x] = bjetmap(num)
if nargin < 1,num = size(colormap,1); end
pt = [0.00 0.00 0.00 0.50;...
      0.13 0.00 0.10 1.00;...
      0.37 0.00 1.00 1.00;...
      0.50 1.00 1.00 1.00;...
      0.57 1.00 1.00 0.00;...
      0.88 1.00 0.00 0.00;...
      1.00 0.50 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = bjetmapinv(num) %#ok<*DEFNU>
if nargin < 1,num = size(colormap,1); end
x = bjetmap(num);
x = flipud(x);
return


function [x] = bjetmapx(num)
if nargin < 1,num = size(colormap,1); end
pt = [0.00 0.00 0.80 0.00;...
      0.15 0.00 0.30 0.00;...
      0.20 0.00 0.00 0.50;...
      0.27 0.00 0.10 1.00;...
      0.44 0.00 1.00 1.00;...
      0.50 1.00 1.00 1.00;...
      0.54 1.00 1.00 0.00;...
      0.73 1.00 0.00 0.00;...
      0.80 0.50 0.00 0.00;...
      0.85 0.50 0.00 0.50;...
      1.00 0.90 0.00 0.85];
x = feval('fleximap',num,pt);
return


function [x] = bjetmapxinv(num)
if nargin < 1,num = size(colormap,1); end
x = bjetmapx(num);
x = flipud(x);
return 


function [x] = carbonemap(num)
if ~exist('num','var'), num = 16; end
x = carbmap(num);
return


function [x] = carbmap(num)
if nargin<1,num = size(colormap,1); end
if num<16, fprintf('Cannot generate map with %d elements\n',num); x = []; return; end
num = num-1;
pt = [0.000                  mybgcolor(); ...
      1/num                  0.10 0.00 0.30; ...  % dark blue
      0.125                  0.33 0.06 0.73; ...  % blue-purple
      round(0.25*num-1)/num  0.50 0.45 1.00; ...  % light blue
      round(0.25*num)/num    0.00 0.40 0.00; ...  % green
      0.325                  0.00 0.65 0.00; ...  % mid-green
      round(0.475*num-1)/num 0.70 0.90 0.70; ...  % light green
      round(0.475*num)/num   0.70 0.70 0.70; ...  % gray
      0.575                  1.00 1.00 0.00; ...  % light yellow
      0.675                  0.95 0.60 0.10; ...  % yellowish orage
      round(0.825*num)/num   0.50 0.30 0.25; ...  % brown
      round(0.825*num+1)/num 1.00 0.30 0.51; ...  % 
      1.000                  0.60 0.00 0.05];
x = feval('fleximap',num+1,pt);
return


function [x] = czmap(num)
if nargin<1,num = 16; end
if num<12, fprintf('Cannot generate map with %d elements\n',num); num = 12; end
num = num-1;
% pt = [0.000                  mybgcolor(); ...
%       1/num                  0.10 1.00 1.00; ...  % cyan blue
%       0.118                  0.00 0.40 1.00; ...  % ocean blue
%       round(0.261*num-1)/num 0.00 0.00 0.75; ...  % dark blue
%       round(0.261*num)/num   0.20 1.00 0.40; ...  % light green
%       0.332                  0.00 0.80 0.00; ...  % green
%       round(0.451*num-1)/num 0.00 0.50 0.00; ...  % dark green
%       round(0.451*num)/num   1.00 1.00 0.20; ...  % light yellow
%       0.5462                 1.00 0.60 0.00; ...  % orange
%       0.6418                 1.00 0.40 0.40; ...  % pink
%       round(0.737*num)/num   0.95 0.00 0.00; ...  % torch read
%       round(0.785*num)/num   0.67 0.00 0.00; ...  % dark red
%       round(0.785*num+1)/num 1.00 0.00 1.00; ...  % magenta
%       (num-1)/num            0.40 0.10 0.60; ...  % dark purple
%       1.000                  1.00 1.00 1.00];
c = (num-1)/num;
d = c/(5-7/num);
m = num-1;
pt = [0.000                          mybgcolor(); ...
      1/num                          0.10 1.00 1.00; ...  % cyan blue
      1/num+(round(d*m)-1)/num       0.00 0.00 0.70; ...  % dark blue
      1/num+(round(d*m))/num         0.00 1.00 0.00; ...  % light green
      1/num+(round(2*d*m)-1)/num     0.00 0.40 0.00; ...  % dark green
      1/num+(round(2*d*m))/num       1.00 1.00 0.00; ...  % yellow
      1/num+(round(3*d*m)-1)/num     1.00 0.50 0.00; ...  % orange
      1/num+(round(3*d*m))/num       1.00 0.00 0.00; ...  % torch red
      1/num+(round(4*d*m)-1)/num     0.40 0.00 0.00; ...  % dark red
      1/num+(round(4*d*m))/num       1.00 0.00 1.00; ...  % magenta
      (num-1)/num                    0.50 0.10 0.60; ...  % purple
      1.000                          1.00 1.00 1.00];
x = feval('fleximap',num+1,pt);
return


function [x] = czmapx(num)
if nargin<1,num = 23; end
if num<16, fprintf('Cannot generate map with %d elements\n',num); num = 16; end
num = num-1;
c = (num-1)/num;
%c = (num-4)/num;
d = c/(7-4/num);
m = num-1;
%m = num-4;
pt = [0.000                          mybgcolor(); ...
      1/num                          0.80 1.00 1.00; ...
      2/num                          0.80 0.60 0.80; ...  % light purple
      2/num+(round(d*m)-1)/num       0.40 0.20 0.40; ...  % dark purple
      2/num+(round(d*m))/num         0.80 0.80 0.60; ...  % light dirty
      2/num+(round(2*d*m)-1)/num     0.40 0.40 0.40; ...  % dark gray
      2/num+(round(2*d*m))/num       0.10 1.00 1.00; ...  % cyan blue
      2/num+(round(3*d*m)-1)/num     0.00 0.00 0.85; ...  % dark blue
      2/num+(round(3*d*m))/num       0.00 1.00 0.00; ...  % light green
      2/num+(round(4*d*m)-1)/num     0.00 0.50 0.00; ...  % dark green
      2/num+(round(4*d*m))/num       1.00 1.00 0.00; ...  % yellow
      2/num+(round(5*d*m)-1)/num     1.00 0.50 0.00; ...  % orange
      2/num+(round(5*d*m))/num       1.00 0.00 0.00; ...  % torch red
      2/num+(round(6*d*m)-1)/num     0.40 0.00 0.00; ...  % dark red
      2/num+(round(6*d*m))/num       1.00 0.00 1.00; ...  % magenta
      (num-1)/(num)                  0.45 0.25 0.90; ...  % purple
      1.000                          1.00 1.00 1.00];
%disp((pt+1/num)*(num))
x = feval('fleximap',num+1,pt);
return



function [x] = czmapxinv(num)
if nargin<1,num = 23; end
if num<16, fprintf('Cannot generate map with %d elements\n',num); num = 16; end
num = num-1;
c = (num-1)/num;
%c = (num-4)/num;
d = c/(7-4/num);
m = num-1;
%m = num-4;
pt = [0.000                          mybgcolor(); ...
      1/num                          0.80 1.00 1.00; ...
      2/num                          0.40 0.20 0.40; ...  % dark purple
	  2/num+(round(2*d*m)-1)/num     0.80 0.60 0.80; ...  % light purple
      2/num+(round(2*d*m))/num       0.00 0.00 0.85; ...  % dark blue
      2/num+(round(3*d*m)-1)/num     0.10 1.00 1.00; ...  % cyan blue
      2/num+(round(3*d*m))/num       0.00 0.50 0.00; ...  % dark green
      2/num+(round(4*d*m)-1)/num     0.00 1.00 0.00; ...  % light green
      2/num+(round(4*d*m))/num       1.00 0.50 0.00; ...  % orange
      2/num+(round(5*d*m)-1)/num     1.00 1.00 0.00; ...  % yellow
      2/num+(round(5*d*m))/num       0.50 0.00 0.00; ...  % dark red
      2/num+(round(6*d*m)-1)/num     1.00 0.00 0.00; ...  % torch red
      2/num+(round(6*d*m))/num       0.50 0.00 0.50; ...  % purple
      num/(num+1)                    1.00 0.00 1.00; ...  % magenta
      1.000                          1.00 1.00 1.00];
x = feval('fleximap',num+1,pt);
return



function [x] = czmapx2(s)
% Shades per color for purple, brown, blue, green, yellow & red
% Purple has 2/3 of the others
if nargin<1,s = 3; end
if s<2 || rem(s, 3)~=0
	fprintf('Poor choice for %d shades/element\n',s);
	fprintf('Good choices are: 30, 15, 6, 3\n');
	%s = 3;
end

num = round(6.6667*s)+2;

n = num - 1;
% Boundaries for each color, plus bgcolor & white
pt = zeros(16, 1);
for i = 1:6
	pt(i*2)   = round((i-1)*s+1)/n;
	pt(i*2+1) = round(i*s)/n;
end
pt(14) = round(6*s+1)/n;
pt(15) = round(6*s+0.66667*s)/n;
pt(16) = 1.0;
%disp((pt+1/n)*(n))	

pt = cat(2, pt, ...
	[mybgcolor(); ...
     0.80 0.60 0.80; ...  % light purple
     0.40 0.20 0.40; ...  % dark purple
     0.80 0.80 0.60; ...  % light dirty
     0.40 0.40 0.40; ...  % dark gray
     0.10 1.00 1.00; ...  % cyan blue
     0.00 0.00 0.90; ...  % dark blue
     0.00 1.00 0.00; ...  % light green
     0.00 0.50 0.00; ...  % dark green
     1.00 1.00 0.00; ...  % yellow
     1.00 0.50 0.00; ...  % orange
     1.00 0.00 0.00; ...  % torch red
     0.50 0.00 0.00; ...  % dark red
     1.00 0.00 1.00; ...  % magenta
     0.56 0.35 1.00; ...  % purple
     1.00 1.00 1.00]);

 x = feval('fleximap',num,pt);
return


function [x, lim] = bzmap(num)
if ~exist('num', 'var'), num = 3; end
num = min(30, num);
x = feval('czmapx2', num);
x = x(2:end, :);
% Check window color, see if I'm using bright/dark theme
n = num*2;
x(1:n,: ) = feval('fleximap', n, ...
	[0    mybgcolor(); ...
	 1.0  0.50 0.60 0.70]);
dz = 15/num;
lim = [-25 75+dz];

return


function [x] = gomap(num)
if nargin < 1,num = size(colormap,1); end
pt = [0.00 0.25 0.10 0.00;...
      0.20 0.85 0.45 0.05;...
      0.50 1.00 1.00 0.65;...
      0.80 0.00 0.85 0.00;...
      1.00 0.00 0.25 0.00];
x = feval('fleximap',num,pt);
return


function [x] = ogmap(num)
if nargin < 1,num = size(colormap,1); end
pt = [0    0.00 0.25 0.00;...
      0.20 0.00 0.85 0.00;...
      0.50 1.00 1.00 1.00;...
      0.80 0.85 0.45 0.05;...
      1.00 0.25 0.10 0.00];
x = feval('fleximap',num,pt);
return


function [x] = refmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.00 1.0 0.5 0.0;...  % 
      0.45 1.0 1.0 0.0;...  %    orange
      0.50 1.0 1.0 1.0;...  %    white
      0.55 0.0 1.0 0.0;...  %    green
      1.00 0.0 0.5 0.0];    % 
x = feval('fleximap',num,pt);
return 


function [x] = rbmap(num)
if nargin < 1,num = size(colormap,1); end
pt = [0    0.0 0.0 0.5;...
      0.20 0.0 0.0 1.0;...
      0.50 1.0 1.0 1.0;...
      0.80 1.0 0.0 0.0;...
      1.00 0.5 0.0 0.0];
x = feval('fleximap',num,pt);
return


function [x] = rgmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0    0.00 0.20 0.00;...
      0.30 0.00 0.80 0.00;...
      0.50 0.85 0.85 0.85;...
      0.70 0.80 0.00 0.00;...
      1.00 0.20 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = rgmapinv(num)
if nargin<1,num = size(colormap,1); end
x = feval('rgmap',num);
idx1 = floor(num/2); idx2 = num-idx1+1;
x = x([idx1:-1:1 idx1+1:idx2-1 num:-1:idx2],:);
return


function [x] = rgmap2(num)
if nargin<1,num = size(colormap,1); end
c = floor(num/2);
pt = [0             0.00 1.00 0.00;...
      (c-2)/(num-1) 0.00 0.40 0.00;...
      (c-1)/(num-1) 0.22 0.33 0.22;...
      (c)/(num-1)   0.40 0.40 0.40;...
      (c+1)/(num-1) 0.33 0.22 0.22;...
      (c+2)/(num-1) 0.45 0.00 0.00;...
      1.000         1.00 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = rgmapw(num)
if nargin<1,num = size(colormap,1); end
c = floor(num/2);
pt = [0             0.00 1.00 0.00;...
      (c-2)/(num-1) 0.00 0.40 0.00;...
      (c-1)/(num-1) 0.45 0.58 0.45;...
      (c  )/(num-1) 0.58 0.45 0.45;...
      (c+1)/(num-1) 0.45 0.00 0.00;...
      1.000         1.00 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = rgmapwe(num)
if nargin<1,num = size(colormap,1); end
anchors = floor([0 0.25 0.25 0.5 0.5 0.5 0.5 0.75 0.75 1.0]*num) + ...
                [1    0    1   -1  0   1   2    0    1   0];
% disp(anchors)
anchors = (anchors(:)-1) / (num-1);
pt = [anchors [0.00 0.00 1.00; ...
               0.00 1.00 1.00; ...
               0.00 1.00 0.00; ...
               0.00 0.40 0.00; ...
               0.45 0.60 0.45; ...
               0.60 0.40 0.40; ...
               0.45 0.00 0.00; ...
               1.00 0.00 0.00; ...
               1.00 0.40 0.00; ...
               1.00 1.00 0.30; ...
               ]];
x = feval('fleximap',num,pt);
return


function [x] = rgmapn(num)
if nargin<1,num = size(colormap,1); end
c = floor(num/2);
pt = [0             0.00 1.00 0.00;...
      (c-2)/(num-1) 0.00 0.45 0.00;...
      (c-1)/(num-1) 0.15 0.27 0.15;...
      (c)/(num-1)   0.27 0.15 0.15;...
      (c+1)/(num-1) 0.45 0.00 0.00;...
      1.000         1.00 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = t7map(num)
if nargin<1,num = size(colormap,1); end
pt = [0     0.00 0.20 0.30;...
      0.350 0.00 0.47 0.35;...
      0.699 0.00 0.75 0.20;...
      0.700 1.00 0.25 1.00;...
      0.950 1.00 1.00 0.00;...
      1.000 1.00 1.00 1.00];
x = feval('fleximap',num,pt);
return


function [x] = t6map(num)
if nargin<1,num = size(colormap,1); end
pt = [0     0.00 0.10 0.25;...
      0.300 0.00 0.32 0.38;...
      0.599 0.00 0.65 0.20;...
      0.600 1.00 0.30 0.90;...
      0.920 1.00 1.00 0.00;...
      1.000 1.00 1.00 1.00];
x = feval('fleximap',num,pt);
return


function [x] = t5map(num)
if nargin<1,num = size(colormap,1); end
pt = [0     0.00 0.10 0.25;...
      0.250 0.00 0.32 0.38;...
      0.499 0.00 0.65 0.20;...
      0.500 1.00 0.20 1.00;...
      0.850 1.00 1.00 0.00;...
      1.000 1.00 1.00 1.00];
x = feval('fleximap',num,pt);
return


function [x] = bgraymap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.000 0.00 0.90 0.10
      0.499 1.00 1.00 1.00;...
      0.500 0.50 0.50 0.50;...
      0.501 0.00 0.00 0.00;...
      1.000 1.00 0.10 0.10];
x = feval('fleximap',num,pt);
return


function [x] = wmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.00 1.0 1.0 1.0;...  % 
      0.25 0.0 1.0 1.0;...  %    cyan
      0.80 0.0 0.0 1.0;...  %    blue
      1.00 0.0 0.0 0.5];    % 
x = feval('fleximap',num,pt);
return 


function [x] = tempmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.000 0.0 0.0 0.8;...  %    blue
	  0.100 0.2 0.5 1.0;...  %    mid blue
      0.200 0.0 1.0 1.0;...  %    fluoroscent cyan
      0.333 0.0 0.7 0.1;...  %    forest green
      0.467 0.5 1.0 0.0;...  %    lemon green
      0.533 1.0 1.0 0.0;...  %    pure yellow
      0.700 1.0 0.25 0.0;...  %    red-orange
      0.733 1.0 0.0 0.0;...  %    husker red
      0.900 0.3 0.0 0.0;...  %    devil's red
      1.000 0.5 0.0 0.6];    %    
x = feval('fleximap',num,pt);
return 


function [x] = tempmapinv(num)
if nargin<1,num = size(colormap,1); end
x = tempmap(num);
x = flipud(x);
return 


function [x] = brmap(num)
if nargin<1,num = size(colormap,1); end
midcolor = 0.7*[1 1 1];
c = 2;
pt = [0.00 0.0 0.0 0.4;...  % Dark blue
      (0.5*num-c-0.01)/num 0.2 0.5 1.0;...  % Light blue
      (0.5*num-c)/num midcolor;...  % Gray
      (0.5*num+c)/num midcolor;...  % Gray
      (0.5*num+c+0.01)/num 1.0 0.4 0.4;...  % Light red
      1.00 0.5 0.0 0.0];    % Dark red
x = feval('fleximap',num,pt);
return


function [x] = rcmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0.00 0.7 0.0 0.3;...  % 
      0.25 1.0 0.5 0.0;...  %    orange
      0.45 1.0 1.0 0.0;...  %    yellow
      0.50 1.0 1.0 1.0;...  %    white
      0.58 0.0 1.0 0.0;...  %    green
      0.78 0.0 0.5 0.3;...  %    dark green
      1.00 0.0 1.0 1.0];    %    cyan
x = feval('fleximap',num,pt);
return 


function [x] = drgmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0    0.00 0.35 0.00;...
      0.15 0.00 0.80 0.00;...
      0.35 0.35 1.00 0.00;...
      0.50 1.00 1.00 1.00;...
      0.65 1.00 0.35 0.35;...
      0.85 1.00 0.00 0.00;...
      1.00 0.45 0.00 0.00];
x = feval('fleximap',num,pt);
return


function [x] = dgbmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0    0.00 0.00 0.65;...
      0.15 0.00 0.00 1.00;...
      0.35 0.00 0.40 1.00;...
      0.50 1.00 1.00 1.00;...
      0.65 0.40 1.00 0.00;...
      0.85 0.00 0.80 0.00;...
      1.00 0.00 0.35 0.00];
x = feval('fleximap',num,pt);
return


function [x] = dbrmap(num)
if nargin<1,num = size(colormap,1); end
pt = [0    0.45 0.00 0.00;...
      0.15 1.00 0.00 0.00;...
      0.35 1.00 0.40 0.40;...
      0.50 1.00 1.00 1.00;...
      0.65 0.00 0.40 1.00;...
      0.85 0.00 0.20 1.00;...
      1.00 0.00 0.00 0.60];
x = feval('fleximap',num,pt);
return


function [x] = drgmapinv(num)
if nargin<1,num = size(colormap,1); end
x = flipud(drgmap(num));
return


function [x] = dgbmapinv(num)
if nargin<1,num = size(colormap,1); end
x = flipud(dgbmap(num));
return


function [x] = dbrmapinv(num)
if nargin<1,num = size(colormap,1); end
x = flipud(dbrmap(num));
return


function [x] = classmap(num)
if nargin<1,num = size(colormap,1); end
% pt = [0    0.00 0.00 0.65;...
%       0.15 0.00 0.00 1.00;...
%       0.30 0.00 0.65 1.00;...
%       0.50 1.00 1.00 1.00;...
%       0.70 0.40 1.00 0.00;...
%       0.85 0.00 0.80 0.00;...
%       1.00 0.00 0.35 0.00];
pt = zeros(6,4);
pt(1:6, 1) = linspace(0, 1, 6);
pt(1,2:4) = [0.00 0.35 0.00];
pt(2,2:4) = [0.80 0.30 0.80];
pt(3,2:4) = [0.80 0.50 0.30];
pt(4,2:4) = [0.40 0.50 1.00];
pt(5,2:4) = [0.50 0.70 1.00];
pt(6,2:4) = [1.00 1.00 1.00];

x = feval('fleximap',num,pt);
return


function [x] = dmap(num)
    if nargin<1, num = 51; end
    den = num - 1;
    pt = [0.00     0.30 0.45 0.50; ...
            9/den  0.60 0.90 1.00; ...
           10/den  0.45 0.20 0.80; ...
           39/den  0.70 0.40 1.00; ...
           40/den  0.50 0.20 0.35; ...
           69/den  1.00 0.50 0.85; ...
           70/den  0.70 0.50 0.15; ...
           99/den  1.00 1.00 0.85; ...
          100/den  1.00 1.00 1.00; ... % 0dB
          129/den  0.00 0.35 1.00; ...
          130/den  0.10 1.00 0.50; ... % 3dB
          159/den  0.00 0.50 0.00; ...
          160/den  1.00 1.00 0.00; ... % 6dB
          189/den  1.00 0.50 0.00; ... 
          190/den  1.00 0.00 0.00; ...
          219/den  0.50 0.00 0.00; ...
          220/den  1.00 0.00 1.00; ...
          249/den  0.50 0.00 0.50; ...
          250/den  1.00 1.00 1.00; ...
          1.00     0.60 1.00 1.00];
    x = feval('fleximap',num,pt);
return


function [x] = rainbowmap(num)
    if nargin < 1, num = 256; end
    pt = [0        0.50 0.00 0.00; ...
          1.0/7    1.00 0.00 0.00; ...
          3.0/7      1.00 1.00 0.00; ...
          3.5/7      0.00 1.00 0.00; ...
          4/7      0.00 0.60 0.00; ...
          6/7      0.00 0.00 1.00; ...
          1        0.50 0.00 1.00];
    x = feval('fleximap', num, pt);
return


function [x] = accentmap(num)
    if nargin < 1, num = 256; end
    pt = [0        0.38 0.82 0.50; ...
          1/7      0.76 0.66 0.82; ...
          2/7      1.00 0.75 0.50; ...
          3/7      1.00 1.00 0.50; ...
          4/7      0.25 0.35 0.70; ...
          5/7      1.00 0.00 0.45; ...
          6/7      0.75 0.33 0.00; ...
          1        0.40 0.40 0.40];
    x = feval('fleximap', num, pt);
return


function [x] = rdpumap(num)
    if nargin < 1, num = 256; end
    pt = [0     1.0000 0.9686 0.9529; ...
          0.5   1.0000 0.0000 0.6000; ...
          1.0   0.3255 0.0000 0.4196];
    x = feval('fleximap', num, pt);
return


function [x] = mapinhexq(map)
if nargin<1, map = rgmap(); end
if size(map,2)~=3
	fprintf('Colormap must be in N x 3');
	x = [];
	return;
end
x = reshape(dec2hex(round(255*map'))',[6 size(map,1)])';
idx = 1;
while(idx<=size(map,1))
	if rem(idx,6)==1
		fprintf('\n');
	end
	fprintf('"%s",',x(idx,1:6));
	idx=idx+1;
end
fprintf('\n');
return; 


function [x] = mapinhex(map)
if nargin<1, map = rgmap(); end
if size(map,2)~=3
	fprintf('Colormap must be in N x 3');
	x = [];
	return;
end
x = reshape(dec2hex(round(255*map'))',[6 size(map,1)])';
idx = 1;
while(idx<=size(map,1))
	if rem(idx,6)==1
		fprintf('\n');
	end
	fprintf('0x%s,',x(idx,1:6));
	idx=idx+1;
end
fprintf('\n');
return; 


function [x] = rapmap(npts,mode)
if nargin < 1,npts = size(colormap,1); mode='lin';
elseif nargin < 2,mode='lin'; end
cromIQ = [0.00, 0.00; ...
         -0.07, 0.13; ...
         -0.11, 0.20; ...
         -0.03, 0.21; ...
          0.23, 0.19; ...
          0.39, 0.11; ...
          0.41,-0.01; ...
          0.39,-0.12; ...
          0.32,-0.19; ...
          0.20,-0.15; ...
          0.00, 0.00];
oldx = 1:11;
newx = linspace(1,11,npts).';
switch (mode)
    case 'lin'
        yiqmap(1:npts,1) = linspace(0,1,npts).';
    case 'exu'
        tmp = linspace(0,1,npts).';
        yiqmap(1:npts,1) = 0.45*tmp+0.55*(exp(1)-exp(1-tmp))/(exp(1)-1);
    case 'exd'
        tmp = linspace(0,1,npts).';
        yiqmap(1:npts,1) = 0.45*tmp+0.55*(exp(tmp)-1)/(exp(1)-1);
    case 'rap'
        Raw = [0.00 0.00 0.00; ...
               0.15 0.15 0.50; ...
               0.30 0.15 0.75; ...
               0.60 0.20 0.50; ...
               1.00 0.25 0.15; ...
               0.90 0.50 0.00; ...
               0.90 0.75 0.10; ...
               0.90 0.90 0.50; ...
               1.00 1.00 1.00];
        newx = linspace(1,9,npts).';
        oldx = 1:9;
        cromIQ = Raw*[0.596,-0.274,-0.322; ...
                      0.211,-0.523, 0.312].';
        yiqmap(1:npts,1) = interp1(oldx,Raw*[0.299,0.587,0.114].',newx,'cubic');
end
% yiqmap(1:npts,2) = interp1(oldx,cromIQ(:,1),newx,'cubic');
% yiqmap(1:npts,3) = interp1(oldx,cromIQ(:,2),newx,'cubic');
yiqmap(1:npts,2) = interp1(oldx,cromIQ(:,1),newx,'pchip');
yiqmap(1:npts,3) = interp1(oldx,cromIQ(:,2),newx,'pchip');
x=yiqmap*[1.0000, 1.0000, 1.0000; ...
          0.9562,-0.2727,-1.1037; ...
          0.6214,-0.6468, 1.7006];
x(x<0)=0; x(x>1)=1;


function [x] = asymap0(npts,lims)
if nargin < 1,npts = 64; end
if nargin < 2,lims = [-9,+4,-9,+5]; end
if length(lims) < 2,fprintf('LIMS must be at least a vector of length 2.\n'); x = []; return; end
if length(lims) < 4,lims = [lims,-9,+5]; end
[minnum,maxnum,minclip,maxclip] = deal(lims(1),lims(2),lims(3),lims(4));
minnum = max(minnum,minclip); maxnum = min(maxnum,maxclip);
rang = maxnum - minnum;
if minnum > 0,fprintf('minnum should be negative.\n'); x = []; return; end
if minclip > 0,fprintf('minclip should be negative.\n'); x = []; return; end
if maxnum < 0,fprintf('maxnum should be positive.\n'); x = []; return; end
if maxclip < 0,fprintf('maxclip should be positive.\n'); x = []; return; end
map_num = maxclip - minclip;
long_map_len = ceil(npts*map_num/rang);
num_clr = ceil(-minclip/map_num*long_map_len);
end_sec = floor(num_clr/4);
clr1 = linspace(0,0.05,num_clr-3*end_sec); clr2 = linspace(0.05,0.4,2*end_sec+1); clr3 = linspace(0.4,1,end_sec+1);
rcurv = [clr1,clr2(2:end),clr3(2:end)];
clr1 = linspace(0,0.20,end_sec); clr2 = linspace(0.20,1,num_clr-2*end_sec+1); clr3 = ones(1,end_sec);
gcurv = [clr1,clr2(2:end),clr3];
clr1 = linspace(0.35,0.7,end_sec); clr2 = linspace(0.7,1,num_clr-end_sec+1);
bcurv = [clr1,clr2(2:end)];
tmp1 = [rcurv',gcurv',bcurv'];
tmp2 = feval('boonlib','rbmap',ceil(2*maxclip/map_num*long_map_len));
tmp2 = tmp2(ceil(size(tmp2,1)/2)+1:end,:);
asymap = [tmp1; tmp2];
stidx = round(long_map_len*(minnum-minclip)/(maxclip-minclip))+1;
x = asymap(stidx:stidx+npts-1,:);
return


function [x] = asymap(num,lims)
if nargin < 1,num = size(colormap,1); end
if nargin < 2,lims = [-9,+4,-9,+5]; end
if length(lims) < 2,fprintf('LIMS must be at least a vector of length 2.\n'); x = []; return; end
if length(lims) < 4,lims = [lims,lims]; end
if lims(1)<lims(3),lims(3) = lims(1); end
if lims(2)>lims(4),lims(4) = lims(2); end
tmp = -lims(3)/(lims(4)-lims(3));
pt = [0           0.00 0.00 0.45;...
      0.25*tmp    0.00 0.30 0.75;...
      0.75*tmp    0.40 1.00 0.90;...
      tmp         1.00 1.00 1.00;...
      0.6+0.4*tmp 1.00 0.00 0.00;...
      1.00        0.50 0.00 0.00];
x = feval('fleximap',ceil(num*(lims(4)-lims(3))/(lims(2)-lims(1))),pt);
stidx = max(floor(num*(lims(1)-lims(3))/(lims(2)-lims(1)))+1,1);
x = x(stidx:stidx+num-1,:);
return


function [x] = zmap(zlim)
if nargin<1,zlim = [0 75]; end
if numel(zlim)==1, zlim = [0 5*zlim]; end
if zlim(1)<0,zlim(1) = 0; end
if zlim(2)>75,zlim(2) = 75; end
%tmp = linspace(zlim(1),zlim(2),16);
tmp = zlim(1):5:zlim(2);
tmp = floor(tmp/5)+1;
std_map = [mybgcolor(); ...
           0.20 1.00 1.00; ...
           0.20 0.60 1.00; ...
           0.00 0.00 1.00; ...
           0.30 1.00 0.00; ...
           0.10 0.80 0.00; ...
           0.00 0.60 0.00; ...
           1.00 1.00 0.00; ...
           1.00 0.75 0.00; ...
           1.00 0.50 0.00; ...
%           1.00 0.30 0.30; ...
           1.00 0.00 0.00; ...
           0.75 0.00 0.00; ...
           0.50 0.00 0.00; ...
           1.00 0.00 0.80; ...
           0.60 0.30 1.00; ...
           1.00 1.00 1.00];
x = std_map(tmp,:);
return




function [x] = zmap2(zlim)
if nargin<1,zlim = [7.5 72.5]; end
if zlim(1)<7.5,zlim(1) = 7.5; end
if zlim(2)>72.5,zlim(2) = 72.5; end
tmp = zlim(1):5:zlim(2);
tmp = floor(tmp/5);
std_map = [  0    0    0; ...
             5  143  143; ...
           128  225   80; ...
            99  185   63; ...
            72  143   48; ...
            44  104   33; ...
            16   63   15; ...
           241  191   16; ...
           240  127   33; ...
           240   15   33; ...
           143    1    2; ...
           177   32  127; ...
           202   63  161; ...
           255  255  255]/255;
x = std_map(tmp,:);
return


function [x] = zmap3(num)
if nargin<1,num = 26; end
if num<12, fprintf('Cannot generate map with %d elements\n',num); num = 12; end
num = num-1;
pt = [0.00   0.00 0.00 0.00; ...
      0.20   0.25 0.20 0.30; ...
      0.30   0.35 0.30 0.43; ...
      0.45   0.60 0.60 0.60; ...
      0.52   0.00 0.75 0.00; ...
      0.58   0.00 0.40 0.00; ...
      0.65   0.70 0.70 0.00; ...
      0.73   0.70 0.30 0.00; ...
      0.78   1.00 0.00 0.00; ...
      0.83   0.55 0.00 0.00; ...
      0.90   0.70 0.00 0.70; ...
      0.93   0.60 0.20 0.80; ...
      1.00   1.00 1.00 1.00];
x = feval('fleximap',num+1,pt);
return


function [x] = zmapn(zlim)
if nargin<1,zlim = [-5 75]; end
if numel(zlim)==1, zlim = [-5 5*zlim]; end
if zlim(1)<-5,zlim(1) = -5; end
if zlim(2)>75,zlim(2) = 75; end
%tmp = linspace(zlim(1),zlim(2),16);
tmp = zlim(1):5:zlim(2);
tmp = floor((tmp-zlim(1))/5)+1;
std_map = [0.00 0.00 0.00; ...
           0.20 0.30 0.30; ...
           0.20 1.00 1.00; ...
           0.20 0.60 1.00; ...
           0.00 0.00 1.00; ...
           0.30 1.00 0.00; ...
           0.10 0.80 0.00; ...
           0.00 0.60 0.00; ...
           1.00 1.00 0.00; ...
           1.00 0.75 0.00; ...
           1.00 0.50 0.00; ...
           1.00 0.00 0.00; ...
           0.75 0.00 0.00; ...
           0.50 0.00 0.00; ...
           1.00 0.00 0.80; ...
           0.60 0.30 1.00; ...
           1.00 1.00 1.00];
x = std_map(tmp,:);
return


function [x] = zmapnx(varargin)
	tmp = zmapn(varargin{:});
%	x = [0.00 0.00 0.00; ...
	x = [0.03 0.06 0.06; ...
	     0.05 0.08 0.08; ...
	     0.10 0.15 0.15; ...
	     tmp(2:end,:)];
return


function [x] = vrmap(lim)
if nargin<1,lim = [-32 32]; end
if lim(1)<-32,lim(1) = 32; end
if lim(2)>32,lim(2) = 32; end
tmp = lim(1):4.5:lim(2);
tmp = floor((tmp+36.5)/4.5);
std_map = [ 24  24 181; ...
            76  76 255; ...
             1 179 179; ...
            76 255 255; ...
             0 179   1; ...
            76 255  76; ...
           102 102 102; ...
            10  10  10; ...
           179 179 179; ...
           255 255  76; ...
           179 179   0; ...
           255  76  76; ...
           179   0   1; ...
           255  76 255; ...
           179   0 179]/255;
x = std_map(tmp,:);
return


function [cmap] = fleximap(num,pt)
if nargin < 1
    num = 15;
    pt = [0    0.5 0.0 0.0;...
          0.20 1.0 0.0 0.0;...
          0.50 1.0 1.0 1.0;...
          0.80 0.0 0.0 1.0;...
          1.00 0.0 0.0 0.5];
end
if any(pt>1),fprintf('PT Matrix cannot have value > 1.0.\n'); return; end
x = linspace(0,1,num); x = x(:);
cmap = interp1(pt(:,1),pt(:,2:4),x,'linear');
cmap(cmap<0) = 0;
cmap(cmap>1) = 1;
return


% ------------------------------ Drawings ------------------------------

function [x] = draw3dbox(lims)
if nargin<1,lims = [0 1 0 1 0 1]; end
% Vertices of each surface, ordered according to right hand rule
nor_sfc_x = [0 0 1 1]';   nor_sfc_y = [1 1 1 1]';    nor_sfc_z = [0 1 1 0]';
sou_sfc_x = [0 1 1 0]';   sou_sfc_y = [0 0 0 0]';    sou_sfc_z = [0 0 1 1]';
eas_sfc_x = [1 1 1 1]';   eas_sfc_y = [0 1 1 0]';    eas_sfc_z = [0 0 1 1]';
wes_sfc_x = [0 0 0 0]';   wes_sfc_y = [1 0 0 1]';    wes_sfc_z = [0 0 1 1]';
top_sfc_x = [0 1 1 0]';   top_sfc_y = [0 0 1 1]';    top_sfc_z = [1 1 1 1]';
bot_sfc_x = [1 1 0 0]';   bot_sfc_y = [1 0 0 1]';    bot_sfc_z = [0 0 0 0]';
% Composite surface
xx = [nor_sfc_x,sou_sfc_x,eas_sfc_x,wes_sfc_x,top_sfc_x,bot_sfc_x];
yy = [nor_sfc_y,sou_sfc_y,eas_sfc_y,wes_sfc_y,top_sfc_y,bot_sfc_y];
zz = [nor_sfc_z,sou_sfc_z,eas_sfc_z,wes_sfc_z,top_sfc_z,bot_sfc_z];
xx = xx*(lims(2)-lims(1))+lims(1);
yy = yy*(lims(4)-lims(3))+lims(3);
zz = zz*(lims(6)-lims(5))+lims(5);
x = patch(xx,yy,zz,'b');
return


function [x] = drawarrow(lims)
if nargin<1,lims = [0 0 0.5]; end
if length(lims)<4,t = pi/4; else, t = lims(4); end
if length(lims)<5,a = 1; else, a = lims(5); end
ha = findobj(gcf,'UserData','ArrowOverlay');
orig_unit = get(gcf,'Unit');
figsize = get(gcf,'Position');
if isempty(ha)
    ha = axes('Unit','Pixel','Position',[1 1 figsize(3:4)]);
elseif length(ha)>1
    delete(ha)
    ha = axes('Unit','Pixel','Position',[1 1 figsize(3:4)]);
else
    axes(ha); set(gca,'Unit','Pixel'); hold on %#ok<*MAXES>
end
axsize = get(gca,'Position');
axaspect = axsize(3)/axsize(4);
a = a/50;
xx = [0 0 lims(3)-a lims(3)-a lims(3) lims(3)-a lims(3)-a 0 0];
yy = [1 -1 -1 -2 0 2 1 1 -1]*0.3*a;
xxp = xx*cos(t)+yy*sin(t)+lims(1);
yyp = (xx*sin(t)-yy*cos(t))*axaspect+lims(2);
tmp = fill(xxp,yyp,'g');
set(ha,'Visible','Off','Unit','Pixel',...
       'XLim',[-0.01 1.01],'XLimMode','Manual',...
       'YLim',[-0.01 1.01],'YLimMode','Manual','UserData','ArrowOverlay');
hold off
set(gcf,'Unit',orig_unit);
x = [ha tmp];
return


function [x] = cspace(cmap)
clf reset
bsizewin(gcf,[1000 500])
len=size(cmap,1);
yiqmap = cmap*[0.299, 0.587, 0.114; ...
               0.596,-0.274,-0.322; ...
               0.211,-0.523, 0.312].';
% Normalize it into [0,1] range
yiqmap(:,2) = (yiqmap(:,2)+0.596)/1.192;
yiqmap(:,3) = (yiqmap(:,3)+0.523)/1.046;

hsvmap = rgb2hsv(cmap);

subplot(2,1,1)
plot(linspace(0,1,len),':','Color',[1 1 1]*0.7);
hold on
x(1)=plot(cmap(:,1),'Color',[1 0 0]);
x(2)=plot(cmap(:,2),'Color',[0 1 0]);
x(3)=plot(cmap(:,3),'Color',[0 0 1]);
x(4)=plot(yiqmap(:,1),'Color',[0.5 0.5 0.5]);
x(5)=plot(yiqmap(:,2),'Color',[0.8 0.8 0.0]);
x(6)=plot(yiqmap(:,3),'Color',[0.4 0.7 1.0]);
x(7)=plot(hsvmap(:,1),'Color',[0.7 0.0 1.0]);
% [legh objh outh outm] = legend(x,'Red Channel','Green Channel','Blue Channel',...
%          'Luminance Y','Chrominance I','Chrominance Q',-1);
[~, objh] = legend(x,'Red Channel','Green Channel','Blue Channel',...
	'Luminance Y','Chrominance I','Chrominance Q','Hue','Location','NorthEastOutside');
tmp = findobj(objh,'Type','Text');
set(tmp,'Color',get(0,'defaultTextColor'));
%tmp = get(0,'defaultAxesXColor');
%set(legh,'XColor',tmp,'YColor',tmp,'ZColor',tmp);
hold off
axis([0.5 len+0.5 0 1]);
xlabel('Colormap Index'); ylabel('Value');
%pos = get(gca,'Position');
pos = [0.07 0.58 0.76 0.34];
set(gca,'Position',[0.07 0.58 0.76 0.34])
if len<50,set(x,'Marker','.'); end

subplot(2,1,2)
img=zeros(8,len,3);
img(1,1:len,1:3)=ind2rgb(1:len,cmap);
img(2,1:len,1:3)=ind2rgb(round(hsvmap(:,1).'*255)+1,hsv(256));
img(3,1:len,1:3)=ind2rgb(round(yiqmap(:,1).'*255)+1,gray(256));
img(4,1:len,1:3)=ind2rgb(round(yiqmap(:,2).'*255)+1,linspace(0,1,256).'*[1 1 0]);
img(5,1:len,1:3)=ind2rgb(round(yiqmap(:,3).'*255)+1,linspace(0,1,256).'*[0.2 0.8 1]);
img(6,1:len,1:3)=ind2rgb(round(cmap(:,1).'*255)+1,[linspace(0,1,256).',zeros(256,2)]);
img(7,1:len,1:3)=ind2rgb(round(cmap(:,2).'*255)+1,[zeros(256,1),linspace(0,1,256).',zeros(256,1)]);
img(8,1:len,1:3)=ind2rgb(round(cmap(:,3).'*255)+1,[zeros(256,2),linspace(0,1,256).']);
image(img)
tmppos=get(gca,'Position');
set(gca,'Position',[pos(1) tmppos(2) pos(3:4)],'YAxisLocation','Right','YTick',1:8,...
    'YTickLabel',{'Colormap','Hue', 'Luminance Y','Chrominance I','Chrominance Q',...
                   'Red Channel','Green Channel','Blue Channel'})
xlabel('Colormap Index');
return


function [xx, yy] = makecircle(params,npts)
if nargin<1,params = [0 0 1]; end
if nargin<2,npts = 90; end
if size(params,2)~=3
    fprintf('Parameter must be in Nx3,e.g. [0 0 1; 0 0 2]\n')
    return
end
x = params(:,1); y = params(:,2); radii = params(:,3);
ang_rad = linspace(0,2*pi,npts);
xx = radii*cos(ang_rad)+repmat(x,[1 npts]);
yy = radii*sin(ang_rad)+repmat(y,[1 npts]);
xx(:,npts+1) = nan;
yy(:,npts+1) = nan;
xx = xx.'; xx = xx(:);
yy = yy.'; yy = yy(:);
return


function [hcir] = circle(params,npts)
if nargin<1,params = [0 0 1]; end
if nargin<2,npts = 90; end
[xx, yy] = makecircle(params,npts);
hold on
hcir = plot(xx,yy,'Color','k');
hold off
return


function [h] = spidergrid(a,sep,nsec)
if nargin<1 || isempty(a), a = gca; end
m = [get(a(1), 'XLim') get(a(1), 'YLim')]; m = max(m(:));
if nargin<2, sep = m/4; end
m = m+2*sep;
if nargin<3, nsec = 12; end
phi = (0:nsec-1)/nsec*2*pi;
b = sin(phi); b = round(b*1e4)*1e-4; % avoid some small misalignment
c = cos(phi); c = round(c*1e4)*1e-4;
xx = kron(b,[nan 0 1])'*m;
xy = kron(c,[nan 0 1])'*m;
r = (sep:sep:m).';
n = length(r);
if (n>0)
	[cx, cy] = makecircle([zeros(n,2) r]);
	if isempty(xx)
		xx = cx;
		yy = cy;
	else
		xx = [cx; xx];
		yy = [cy; xy];
	end
end
h = zeros(length(a),1);
for idx = 1:length(a)
	c = get(a(idx),'Color');
	hold(a(idx),'on')
	h(idx) = plot(a(idx),xx,yy,'Color',1-c);
	hold(a(idx),'off')
end
return;


% ------------------------------ Don't ask ------------------------------
function grayme(h,blend)
if nargin<1, h = gcf; end
if (nargin==1)&&~ishandle(h), blend = h; h = gcf; end
if ~exist('blend','var'), blend = 1; end
if ~ishandle(h), fprintf('Invalid handle object.\n'); return; end
feval('tintme',h,[0.5 0.5],blend);
return


function greenme(h,blend)
if nargin<1, h = gcf; end
if (nargin==1)&&~ishandle(h), blend = h; h = gcf; end
if ~exist('blend','var'), blend = 1; end
if ~ishandle(h), fprintf('Invalid handle object.\n'); return; end
feval('tintme',h,[0.41 0.37],blend);
return


function sepiame(h,blend)
if nargin<1, h = gcf; end
if (nargin==1)&&~ishandle(h), blend = h; h = gcf; end
if ~exist('blend','var'), blend = 1; end
if ~ishandle(h), fprintf('Invalid handle object.\n'); return; end
feval('tintme',h,[0.43 0.53],blend);
return


function tintme(h,adj,blend)
if nargin<1, h = gcf; end
if nargin<2||length(adj)~=2, adj = [0.6 0.4]; end
if ~exist('blend','var'), blend = 1; end
if ~ishandle(h), fprintf('Invalid handle object.\n'); return; end
cmap = get(h,'Colormap');
tmp = rgb2ntsc(cmap);
cmap = rgb2ycbcr(cmap);
cmap(:,1) = tmp(:,1);
cmap(:,2) = blend*adj(1)+(1-blend)*cmap(:,2);
cmap(:,3) = blend*adj(2)+(1-blend)*cmap(:,3);
cmap = ycbcr2rgb(cmap);
clr = get(h,'Color');
clr = rgb2ycbcr(clr); clr(:,2) = 0.4+0.2*adj(1); clr(:,3) = 0.4+0.2*adj(2); clr = ycbcr2rgb(clr);
set(h,'Colormap',cmap,'Color',clr);
%tags = {'Color','BackgroundColor','ForegroundColor','MarkerFaceColor','MarkerEdgeColor','XColor','YColor'};
tags = {'Color','MarkerFaceColor','MarkerEdgeColor','XColor','YColor','FaceColor'};
ha = get(h,'Children');
for idx=1:length(ha)
    for itag=1:length(tags)
        if isprop(ha(idx),tags{itag})
            clr = get(ha(idx),tags{itag});
            if isnumeric(clr)
                %fprintf('Prop = %s  Val = %s\n',tags{itag},num2str(tmp));
				clr = rgb2ycbcr(clr);
				clr(:,2) = blend*adj(1)+(1-blend)*clr(:,2);
				clr(:,3) = blend*adj(2)+(1-blend)*clr(:,3);
				clr = ycbcr2rgb(clr);
                set(ha(idx),tags{itag},clr);
            end
        end
    end
    % Second Level
    hb = get(ha(idx),'Children');
	if strcmp(get(ha(idx),'Type'),'axes')
    	hb = [hb; get(ha(idx),'XLabel');
        get(ha(idx),'YLabel');
        get(ha(idx),'Title')]; %#ok<AGROW>
	end
    for jdx=1:length(hb)
        for itag=1:length(tags)
            if isprop(hb(jdx),tags{itag})
                clr = get(hb(jdx),tags{itag});
                if isnumeric(clr)
                    %fprintf('Prop = %s  Val = %s\n',tags{itag},num2str(tmp));
					clr = rgb2ycbcr(clr);
					clr(:,2) = blend*adj(1)+(1-blend)*clr(:,2);
					clr(:,3) = blend*adj(2)+(1-blend)*clr(:,3);
					clr = ycbcr2rgb(clr);
                    set(hb(jdx),tags{itag},clr);
                end
            end
        end
    end
end
% Special image type
ha = findobj(h,'Type','Image');
for idx=1:length(ha)
    tmp = get(ha(idx),'CData');
    if size(tmp,3)==3
        % The RGB type
		clr = rgb2ycbcr(tmp);
		clr(:,2) = blend*adj(1)+(1-blend)*clr(:,2);
		clr(:,3) = blend*adj(2)+(1-blend)*clr(:,3);
		clr = ycbcr2rgb(clr);
        set(ha(idx),'CData',clr);
    end
end
return


function [fname] = choosefile(arg1,arg2,arg3)
if nargin<2,arg2 = '*'; end
if nargin<1,arg1 = './'; end
if iscell(arg1), path_choice = arg1;
elseif exist(arg1,'dir')==7, path_choice = {arg1};
elseif exist(arg1,'file')==2, path_choice = feval(arg1);
elseif ~exist(arg1,'dir') && ~exist(arg1,'file')
    fprintf('Path ''%s'' does not exist.\n', arg1);
    fname = [];
    return;
end
if exist('arg3','var') && numel(arg3)~=2, clear arg3; end
% Initialize an empty fname for early returns
fname = [];
if length(path_choice)>1
	if exist('arg3','var')&&((arg3(1)>0)&&(arg3(1)<=length(path_choice)))
		tmp = arg3(1);
		% fprintf('Pre-selected Folder #%d: %s\n',tmp,char(path_choice(tmp)));
	else
		for idx = 1:length(path_choice)
			fprintf(' % 3d. %s\n',idx,char(path_choice{idx}));
		end
		fprintf('\n');
		tmp = input(['Selection (1-',num2str(idx),') : ']);
		if isempty(tmp) || ~ismember(tmp,1:idx)
			fprintf('No valid selection,action cancelled !\n\n');
			return
		end
	end
	dirname = char(path_choice(tmp));
else
	% fprintf('Only one folder.  Automatically selected.\n');
    dirname = char(path_choice);
end
if ~exist(dirname,'dir')
    fprintf('Directory does not exists, double check the path(s).\n');
    fname = [];
    return;
end
if dirname(end)~='/', dirname = [dirname,'/']; end
flist = dir([dirname,arg2]);
if isempty(flist)
    fprintf('Folder ''%s'' is empty.\n', dirname);
    return;
end
% Take out . and ..
names = {flist.name}';
keep_idx = ~ismember(names,{'.','..','.DS_Store'});
flist = flist(keep_idx);
if isempty(flist)
	fprintf('No file matches the criteria.\n');
	fname = [];
	return;
end

if ~exist('arg3','var')
	flist = feval('showflist',flist);
else
	flist = feval('showflist',flist,1);
end

% Selecting the file
if exist('arg3','var')
	if ((arg3(2)>0)&&(arg3(2)<=length(flist)))
		tmp = arg3(2);
	elseif (arg3(2)<=0)
		tmp = max(length(flist)+arg3(2),1);
	else
		%fprintf('%s: Empty filename.\n',mfilename)
		fname = [];
		return;
	end
	fprintf('Pre-selected File #%d: %s\n',tmp,char(flist(tmp)));
else
	idx = length(flist);
	tmp = input(['\nSelection (1-',num2str(idx),') : ']);
	if isempty(tmp) || ~ismember(tmp,1:idx)
		fprintf('No valid selection,action cancelled !\n\n');
		fname = [];
		return
	end
	fprintf('\n');
end

% Construct the full path filename
if strcmp(dirname,'.')
    fname = char(flist(tmp));
else
    fname = [dirname,char(flist(tmp))];
end
return


function x = showflist(flist,quiet)
if isempty(flist), x = []; return; end
if ~exist('quiet','var'), quiet = 0; end
if ~isfield(flist,'bytes')&&~isfield(flist,'name')
    fprintf('Input must be something returned by dir().\n')
    x = 0;
    return
end
% Collect files, display them.
fbyte = [flist.bytes]';
flist = {flist.name}';

tmp = regexp(flist,'(?<=\D*)\d+','match');
if ~isempty(tmp)
	% tmp = cat(1,tmp{:});   % <-- Array with numbers for each entry
	% Convert them into numbers
	tmp2 = zeros(size(tmp));
	for idx=1:numel(tmp2)
		tmp3 = str2double(char(tmp{idx}));
		tmp2(idx,1:length(tmp3)) = tmp3;
	end
	% Try to sort them in a smart way, using those numbers
	[~, tmp] = sortrows(tmp2);
	flist = flist(tmp);
% Try to find indexing scheme, this is just special case for my xpol project
% if ~isempty(strfind(flist{1},'_P'))|~isempty(strfind(flist{1},'_H'))
% 	dd = zeros(size(fbyte,1),1);
% 	for idx = 1:length(dd)
% 		ii = strfind(flist{idx},'_P');
% 		if isempty(ii), ii = strfind(flist{idx},'_H'); end
% 		tmp = flist{idx}(ii+2:length(flist{idx})-4);
% 		keyboard
% 		dd(idx) = str2num(tmp);
% 	end
% 	[dds tmp] = sort(dd);
% 	flist = flist(tmp);
else
 	[flist, tmp] = sortrows(flist);
end
if ~quiet
	fbyte = fbyte(tmp);
	if length(fbyte)<1,msize = 0; else, msize = mean(fbyte); end
	if msize<1e3,funits = 'B';
	elseif msize>1e3&&msize<1e6, fbyte = fbyte/1024; funits = 'KB';
	elseif msize>1e6&&msize<1e9, fbyte = fbyte/1024/1024; funits = 'MB';
	elseif msize>1e9&&msize<1e12, fbyte = fbyte/1024/1024/1024; funits = 'GB';
    else, fprintf('Filesize is rediculously big at the year of 2005.\n'); x = []; return;
	end
	tmp = size(char(flist));
	msize = ceil(log10(max(fbyte)+1));
	ndig = ceil(log10(length(flist)+1));
	% Compact template
	ftemplate = ['%',num2str(ndig),'d. %',num2str(tmp(2)),'s %',num2str(msize+3),'.2f ',funits,'\n'];
	if ((tmp(2)+ndig+msize+3+5)<39)
		% Assume 79 screen width
		% ndig digits + 4 spaces; 
		% tmp(2) = filename's length; msize = filesize digits
		leftover = 79-2*(ndig+4+tmp(2)+msize+3+length(funits));
		pad = char(repmat(' ',[1 max(floor(leftover/3),1)]));
		ftemplate = [pad,ftemplate];
	else
		% Can't squeeze the info in two columns ==> space them out a bit
		ftemplate = ['%',num2str(ndig),'d.  %',num2str(tmp(2)),'s   %',num2str(msize+3),'.2f ',funits,'\n'];
		leftover = 79-(ndig+7+tmp(2)+msize+3+length(funits));
		pad = char(repmat(' ',[1 floor(leftover/2)]));
		ftemplate = [pad,ftemplate];
	end
	if length(flist)==1
		idx = 1;
		fprintf(ftemplate,idx,char(flist(idx)),fbyte(idx)) %#ok<PRTCAL>
	else
		if ((tmp(2)+ndig+msize+3+5)<39)
			half_idx = round(length(flist)/2);
			for idx=1:half_idx-1+double(mod(length(flist),2)==0)
				fprintf(ftemplate(1:end-2),idx,char(flist(idx)),fbyte(idx)); % Left column
				fprintf(ftemplate,half_idx+idx,char(flist(half_idx+idx)),fbyte(half_idx+idx)); % Right column
			end
			if mod(length(flist),2)~=0
				fprintf(ftemplate,half_idx,char(flist(half_idx)),fbyte(half_idx));
			end
		else
			for idx=1:length(flist)
				fprintf(ftemplate,idx,char(flist(idx)),fbyte(idx));
			end
		end
	end
end
x = flist;
return


function [yesno] = aintersect(fh)
yesno = false;
if length(fh)<2,return; end
if nargin<1
    pos_set = get(get(0,'Children'),'Position');
    pos_set = reshape([pos_set{:}],4,length(pos_set)).';
else
    pos_set = get(fh,'Position');
    pos_set = reshape([pos_set{:}],4,length(pos_set)).';
end
xbuff = 1; ybuff = 1;
for iset = 1:size(pos_set,1)
    for iobj = 1:size(pos_set,1)
        if iset~=iobj
            tmp = [((pos_set(iobj,1)>=pos_set(iset,1)&...
                     pos_set(iobj,1)<=pos_set(iset,1)+pos_set(iset,3))&...
                    (pos_set(iobj,2)>=pos_set(iset,2)&...
                     pos_set(iobj,2)<=pos_set(iset,2)+pos_set(iset,4))),...
                   ((pos_set(iobj,1)+pos_set(iobj,3)+xbuff>=pos_set(iset,1)&...
                     pos_set(iobj,1)+pos_set(iobj,3)+xbuff<=pos_set(iset,1)+pos_set(iset,3))&...
                    (pos_set(iobj,2)>=pos_set(iset,2)&...
                     pos_set(iobj,2)<=pos_set(iset,2)+pos_set(iset,4))),...
                   ((pos_set(iobj,1)>=pos_set(iset,1)&...
                     pos_set(iobj,1)<=pos_set(iset,1)+pos_set(iset,3))&...
                    (pos_set(iobj,2)+pos_set(iobj,4)+ybuff>=pos_set(iset,2)&...
                     pos_set(iobj,2)+pos_set(iobj,4)+ybuff<=pos_set(iset,2)+pos_set(iset,4))),...  
                   ((pos_set(iobj,1)+pos_set(iobj,3)+xbuff>=pos_set(iset,1)&...
                     pos_set(iobj,1)+pos_set(iobj,3)+xbuff<=pos_set(iset,1)+pos_set(iset,3))&...
                    (pos_set(iobj,2)+pos_set(iobj,4)+ybuff>=pos_set(iset,2)&...
                     pos_set(iobj,2)+pos_set(iobj,4)+ybuff<=pos_set(iset,2)+pos_set(iset,4)))];
            if any(tmp), yesno = true; return; end
%         else
%             tmp = logical([0 0 0 0]);
        end
    end
end
return


function [x] = assignfig(tag,forcenew)
if nargin<2, forcenew=0; end
if nargin<1, fprintf('I cannot comply, no Tag for search.\n'); return; end
if ~ischar(tag), fprintf('Tag must be a string that was assigned to a figure.\n'); return; end
stidx = 3000;
interv = 20;
if (forcenew>1)
	% If a specific number was given, then just use that as figure handle
	x = forcenew;
	return
end
% Find the object with the match tag
tmp = findobj('Tag',tag);
if ~isempty(tmp), tmp = tmp(ishandle(tmp)); end
if ~isempty(tmp)
    % There is at least 1 Figure
	if (forcenew==1)
        % If forcenew, try to open up a new figure that does not overlap
        ifig = round(stidx+interv*rand(1));
        while ismember(ifig,tmp), ifig = ifig+1; end
        x = ifig;
        return
	elseif (forcenew==0)
        % Don't open up new figure, just use to existing one
        x = tmp(1);
        return
	end
else
    % If figure doesn't exist, just give any handle
    x = round(stidx+interv*rand(1));
    return
end


function [yn] = isoverlap(L,R)
% L = (xmin xmax ymin ymax)
if numel(L)~=4
    fprintf('The first input must be a 4-element vector only.\n');
    yn = -1;
    return
end
if size(R,1)~=4
    if numel(R)==4
        R = R(:);
    else
        fprintf('Sorry I can''t compute.\n');
        yn = -1;
        return
    end
end
yn = ( (L(1)>R(1,:))&(L(1)<R(2,:))&...
       (L(3)>R(3,:))&(L(3)<R(4,:)) )|...
     ( (L(1)>R(1,:))&(L(1)<R(2,:))&...
       (L(4)>R(3,:))&(L(4)<R(4,:)) )|...
     ( (L(2)>R(1,:))&(L(2)<R(2,:))&...
       (L(3)>R(3,:))&(L(3)<R(4,:)) )|...
     ( (L(2)>R(1,:))&(L(2)<R(2,:))&...
       (L(4)>R(3,:))&(L(4)<R(4,:)) )|...
     ( (R(1,:)>L(1))&(R(1,:)<L(2))&...
       (R(3,:)>L(3))&(R(3,:)<L(4)) )|...
     ( (R(1,:)>L(1))&(R(1,:)<L(2))&...
       (R(4,:)>L(3))&(R(4,:)<L(4)) )|...
     ( (R(2,:)>L(1))&(R(2,:)<L(2))&...
       (R(3,:)>L(3))&(R(3,:)<L(4)) )|...
     ( (R(2,:)>L(1))&(R(2,:)<L(2))&...
       (R(4,:)>L(3))&(R(4,:)<L(4)) );
return


function singlebar()
pt = cat(2, [0; 1], hsv2rgb([0.6 0.1 0.8; 0.6 0.9 0.8]));
clf
n = 200;
cmap = boonlib('fleximap',n,pt);
colormap(cmap)
imagesc((0:n-1)/(n-1),1,(0:n-1)/(n-1))
set(gca,'position',[0.025 0.33 0.95 0.5],'FontSize',16,'YTick',[],...
    'TickLength',[0 0])
bsizewin(gcf,[1000 80])


function hsvdemo()
tmp = findobj('Tag', 'HSVDemo');
newfig = false;
if isempty(tmp) || ~ishandle(tmp)
    newfig = true;
else
    FIG = get(tmp, 'UserData');
end
if ~exist('FIG','var') || ~all(ishandle(FIG.ha))
    newfig = true;
end
if newfig
	FIG.hf = figure(1001);
	clf
	FIG.n = 200;
	FIG.ii = (0:FIG.n-1)/(FIG.n-1);
	FIG.scl = 3;
	FIG.w = FIG.scl*FIG.n;
	FIG.h = FIG.scl*8;
	FIG.s = 45;
	FIG.ox = 50;
	FIG.oy = 50;
	FIG.hsv = [0.5 0.5 0.5];
	FIG.rgb = [0.5 0.5 0.5];
    FIG.format = ' %.3f ';
	FIG.format2 = ['\\color{magenta}H:',FIG.format,...
		           ' \\color{cyan}S:',FIG.format,...
				   ' \\color{yellow}V:',FIG.format];
	FIG.format3 = ['\\color{red}R:',FIG.format,...
		           ' \\color{green}G:',FIG.format,...
				   ' \\color[rgb]{.1 .4 1}B:',FIG.format];

	FIG.ha = axes('Unit','Pixel','Position',[FIG.ox FIG.oy+2*(FIG.s+FIG.h) FIG.w FIG.h]);
	tmp = ind2rgb(1:FIG.n,hsv(FIG.n));
	FIG.him = image(FIG.ii,0.5,tmp);
	hold on
	FIG.hl = plot(FIG.ha(1),[0.5 0.5],[0 1],'w','LineWidth',2);
	hold off
	FIG.ht = text(0.5,-0.2,sprintf(FIG.format,FIG.hsv(1)));

	FIG.ha(2) = axes('Unit','Pixel','Position',[FIG.ox FIG.oy+FIG.s+FIG.h FIG.w FIG.h]);
	cmap = fleximap(FIG.n,cat(2,[0; 1],hsv2rgb([FIG.hsv(1) 0 FIG.hsv(3); FIG.hsv])));
	tmp = ind2rgb(1:FIG.n,cmap);
	FIG.him(2) = image(FIG.ii,0.5,tmp);
	hold on
	FIG.hl(2) = plot(FIG.ha(2),[0.5 0.5],[0 1],'w','LineWidth',2);
	hold off
	FIG.ht(2) = text(0.5,-0.2,sprintf(FIG.format,FIG.hsv(2)));

	FIG.ha(3) = axes('Unit','Pixel','Position',[FIG.ox FIG.oy FIG.w FIG.h]);
	tmp = ind2rgb(1:FIG.n,gray(FIG.n));
	FIG.him(3) = image(FIG.ii,0.5,tmp);
	hold on
	FIG.hl(3) = plot(FIG.ha(3),[0.5 0.5],[0 1],'w','LineWidth',2);
	hold off
	FIG.ht(3) = text(0.5,-0.2,sprintf(FIG.format,FIG.hsv(3)));

	FIG.hc = axes('Unit','Pixel','Position',...
		[FIG.ox+FIG.w+FIG.s FIG.oy (3*FIG.h+2*FIG.s)*[1 1]]);
	FIG.hcim = image(0,0,ind2rgb(1,hsv2rgb(FIG.hsv)));
	set(FIG.hc,'XTick',[],'YTick',[]);
	FIG.hct = text(0.1,-0.52,sprintf(FIG.format2,FIG.hsv));
	FIG.hct(2) = text(0.1,0.52,sprintf(FIG.format3,FIG.rgb));

	feval('bsizewin',FIG.hf,[FIG.w+3*FIG.h+3*FIG.s+100 3*FIG.h+2*FIG.s+100]);
	FIG.mouse_down = false;
	FIG.bar = 0;
	FIG.code = ['tmp = get(gcf, ''CurrentPoint''); ',...
         'tmp = tmp(1)-FIG.ox; tmp = tmp/FIG.n/FIG.scl; ',...
         'tmp = max(0,min(1,tmp)); ',...
         'FIG.hsv(FIG.bar) = tmp; ',...
         'set(FIG.hl(FIG.bar),''XData'',[tmp tmp]); ',...
         'tmp = cat(2,[0; 1],hsv2rgb([FIG.hsv(1) 0 FIG.hsv(3); FIG.hsv(1) 1 FIG.hsv(3)])); ',...
         'cmap = feval(''boonlib'',''fleximap'',FIG.n,tmp); ',...
         'set(FIG.him(2),''CData'',ind2rgb(1:FIG.n,cmap)); ',...
		 'FIG.rgb = hsv2rgb(FIG.hsv); ',...
		 'set(FIG.hcim,''CData'',ind2rgb(1,FIG.rgb)); ',...
		 'for ii=1:3, ',...
		 '  tmp = get(FIG.ht(ii),''Position''); tmp(1) = FIG.hsv(ii); ',...
		 '  set(FIG.ht(ii),''Position'',tmp,''String'',sprintf(FIG.format,FIG.hsv(ii))); ',...
		 'end; ',...
		 'set(FIG.hct(1),''String'',sprintf(FIG.format2,FIG.hsv)); ',...
		 'set(FIG.hct(2),''String'',sprintf(FIG.format3,FIG.rgb)); ',...
         ];
	set(FIG.ht,'HorizontalAlignment','Center','VerticalAlignment','Bottom','EdgeColor',[1 1 1]);
	set(FIG.hct(1),'HorizontalAlignment','Center','VerticalAlignment','Bottom');
	set(FIG.hct(2),'HorizontalAlignment','Center','VerticalAlignment','Top',...
		'FontWeight','Bold');
    set(FIG.hf,'Menubar','None','Tag', 'HSVDemo','UserData',FIG)

    set(FIG.ha,'XLim',[0 1],'TickLength',[0 0],'YTick',[]);
    if (ispc)
        set(FIG.ha,'FontSize',10);
    else
        set(FIG.ha,'FontSize',14);
    end
    set([FIG.him FIG.hl FIG.ht],'ButtonDownFcn',...
        ['FIG = get(gcf,''UserData''); ',...
         'FIG.mouse_down = true; ',...
         '[tmp FIG.bar] = ismember(gca, FIG.ha); ',...
         'if ~tmp, [~, FIG.bar] = ismember(gco, FIG.hl); end; ',...
         'if ~tmp, [~, FIG.bar] = ismember(gco, FIG.ht); end; ',...
         FIG.code,...
         ]);
    set(FIG.hf, ...
        'WindowButtonUpFcn' ,...
        ['if exist(''FIG'',''var'') && strcmp(get(FIG.hf,''Tag''),''HSVDemo'') && FIG.mouse_down, ',...
		 'FIG.mouse_down = false; ',...
         'FIG.bar = 0; ',...
         'set(FIG.hf,''UserData'',FIG); clear all; ',...
		 'end; '],...
        'WindowButtonMotionFcn',...
        ['if exist(''FIG'',''var'') && FIG.mouse_down, ',...
         FIG.code,...
         'end; ',...
		 ]);
	 clear FIG
end


function mydefault()
set(0,'defaultFigureColor',[0.15 0.15 0.15],'defaultAxesColor',[1 1 1]*0.3,...
      'defaultAxesXColor',[1 1 1],'defaultAxesYColor',[1 1 1],'defaultAxesZColor',[1 1 1],...
      'defaultAxesYDir','Normal',...
	  'defaultAxesTickLength',[0.005 0.005],...
      'defaultTextColor',[1 1 1],...
      'defaultTextFontSize',12,'defaultAxesFontSize',12,...
      'defaultAxesColorOrder',[0.5 1 1 ; 1 0.7 0.0; 0.5 1 0.3; 0 0.8 0; 0 0.3 1; 1 0.6 0; 0 0.8 0.7; 0.6 0.2 1; 0.6 0.6 0.6],...
      'defaultFigureColormap',rapmap(64),'defaultFigureInvertHardCopy','off');

function mydefault2()
set(0,'defaultFigureColor',[1 1 1],'defaultAxesColor',[1 1 1],...
      'defaultAxesXColor',[0 0 0],'defaultAxesYColor',[0 0 0],'defaultAxesZColor',[0 0 0],...
      'defaultAxesYDir','Normal',...
	  'defaultAxesTickLength',[0.005 0.005],...
      'defaultTextColor',[0 0 0],...
      'defaultTextFontSize',12,'defaultAxesFontSize',12,...
      'defaultAxesColorOrder',[0.35 0.50 1.0; 1.0 0.65 0.0; 0.8 0.3 0.8; 0 0.55 0; 0 0.6 0.6; 0.2 0.4 1; 1 0.2 1],...
      'defaultFigureColormap',rapmap(64));
% [1 0 0; 0 0.8 0; 0 0.3 1; 1 0.6 0; 0 0.8 0.7; 0.6 0.2 1; 0.6 0.6 0.6]

function mydefault3()
set(0,'defaultFigureColor',[0.15 0.20 0.20],'defaultAxesColor',[0 0 0],...
      'defaultAxesXColor',[1 1 1],'defaultAxesYColor',[1 1 1],'defaultAxesZColor',[1 1 1],...
      'defaultAxesYDir','Normal',...
	  'defaultAxesTickLength',[0.005 0.005],...
      'defaultTextColor',[1 1 1],...
      'defaultTextFontSize',12,'defaultAxesFontSize',12,...
      'defaultAxesColorOrder',[1 0 0; 0 0.8 0; 0 0.3 1; 1 0.6 0; 0 0.8 0.7; 0.6 0.2 1; 0.6 0.6 0.6],...
      'defaultFigureColormap',rapmap(64),'defaultFigureInvertHardCopy','off');

  
function v = figbrightness()
if ~isempty(gcbf)
	color_hsv = rgb2hsv(get(gcbf, 'Color'));
else
	color_hsv = rgb2hsv(get(0, 'DefaultFigureColor'));
end
v = color_hsv(3);
return
	  

function ind = nwsd2ind(data)
ind = zeros(size(data));
lvl = [-4 -2 -0.5 0 0.25 0.5 1.0 1.5 2 2.5 3 4 5 6 8];
for ii = 1:numel(lvl)
    mask = data >= lvl(ii);
    ind(mask) = ii;
end


function ind = nwsr2ind(data)
ind = zeros(size(data));
lvl = [0.2 0.45 0.65 0.75 0.80 0.85 0.90 0.93 0.95 0.96 0.97 0.98 0.99 1.00];
for ii = 1:numel(lvl)
    mask = data >= lvl(ii);
    ind(mask) = ii;
end


function c = nwsv7map()
c = [hex2dec('00') hex2dec('00') hex2dec('00'); ...
     hex2dec('00') hex2dec('E0') hex2dec('FF'); ...
     hex2dec('00') hex2dec('80') hex2dec('FF'); ...
     hex2dec('32') hex2dec('00') hex2dec('96'); ...
     hex2dec('00') hex2dec('FB') hex2dec('90'); ...
     hex2dec('00') hex2dec('BB') hex2dec('99'); ...
     hex2dec('00') hex2dec('8F') hex2dec('00'); ...
     hex2dec('CD') hex2dec('C9') hex2dec('9F'); ...
     hex2dec('76') hex2dec('76') hex2dec('76'); ...
     hex2dec('F8') hex2dec('87') hex2dec('00'); ...
     hex2dec('FF') hex2dec('CF') hex2dec('00'); ...
     hex2dec('FF') hex2dec('FF') hex2dec('00'); ...
     hex2dec('AE') hex2dec('00') hex2dec('00'); ...
     hex2dec('D0') hex2dec('70') hex2dec('00'); ...
     hex2dec('FF') hex2dec('00') hex2dec('00'); ...
     hex2dec('77') hex2dec('00') hex2dec('7D')] / 255;
 return


function c = nwsvmap()
c = [hex2dec('00') hex2dec('00') hex2dec('00'); ...
     hex2dec('00') hex2dec('E0') hex2dec('FF'); ...
     hex2dec('00') hex2dec('BB') hex2dec('00'); ...
     hex2dec('00') hex2dec('8F') hex2dec('00'); ...
     hex2dec('F8') hex2dec('87') hex2dec('00'); ...
     hex2dec('FF') hex2dec('CF') hex2dec('00'); ...
     hex2dec('FF') hex2dec('00') hex2dec('00'); ...
     hex2dec('77') hex2dec('00') hex2dec('7D')] / 255;
return

     
function c = nwsdmap()
c = [hex2dec('00') hex2dec('00') hex2dec('00'); ...
     hex2dec('40') hex2dec('40') hex2dec('40'); ...
     hex2dec('9C') hex2dec('9C') hex2dec('9C'); ...
     hex2dec('C9') hex2dec('C9') hex2dec('C9'); ...
     hex2dec('8C') hex2dec('78') hex2dec('B4'); ...
     hex2dec('00') hex2dec('00') hex2dec('98'); ...
     hex2dec('23') hex2dec('98') hex2dec('D3'); ...
     hex2dec('44') hex2dec('FF') hex2dec('D2'); ...
     hex2dec('57') hex2dec('DB') hex2dec('56'); ...
     hex2dec('FF') hex2dec('FF') hex2dec('60'); ...
     hex2dec('FF') hex2dec('90') hex2dec('45'); ...
     hex2dec('DA') hex2dec('00') hex2dec('00'); ...
     hex2dec('AE') hex2dec('00') hex2dec('00'); ...
     hex2dec('F7') hex2dec('82') hex2dec('BE'); ...
     hex2dec('FF') hex2dec('FF') hex2dec('FF'); ...
     hex2dec('77') hex2dec('00') hex2dec('7D')] / 255;
return
 

function c = nwsrmap()
c = [hex2dec('00') hex2dec('00') hex2dec('00'); ...
     hex2dec('95') hex2dec('94') hex2dec('9C'); ...
     hex2dec('16') hex2dec('14') hex2dec('8C'); ...
     hex2dec('09') hex2dec('02') hex2dec('D9'); ...
     hex2dec('89') hex2dec('87') hex2dec('D6'); ...
     hex2dec('5C') hex2dec('FF') hex2dec('59'); ...
     hex2dec('8B') hex2dec('CF') hex2dec('02'); ...
     hex2dec('FF') hex2dec('FB') hex2dec('00'); ...
     hex2dec('FF') hex2dec('C4') hex2dec('00'); ...
     hex2dec('FF') hex2dec('89') hex2dec('03'); ...
     hex2dec('FF') hex2dec('2B') hex2dec('00'); ...
     hex2dec('E3') hex2dec('00') hex2dec('00'); ...
     hex2dec('A1') hex2dec('00') hex2dec('00'); ...
     hex2dec('97') hex2dec('05') hex2dec('56'); ...
     hex2dec('FA') hex2dec('AC') hex2dec('D1'); ...
     hex2dec('77') hex2dec('00') hex2dec('7D')] / 255;
return


function h = nwsdbar()
caxis([0 16])
h = colorbar;
set(h, 'YTick', 0:15, 'YTickLabel', ...
    {'< TH', '-4.0', '-2.0', '-0.5', '0.0', '0.25', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '4.0', '5.0', '6.0', 'RF'});
return


function h = nwsrbar()
caxis([0 16])
h = colorbar;
set(h, 'YTick', 0:15, 'YTickLabel', ...
    {'< TH', '0.2', '0.45', '0.65', '0.75', '0.80', '0.85', '0.90', '0.93', '0.95', '0.96', '0.97', '0.98', '0.99', '1.00', 'RF'});
return

function pdata = pdata_setup()
pdata = struct('value', zeros(2, 2), ...
    'title', 'Title', ...
    'colormap', 'zmap', ...
    'clim', [0 16]);
return

% pdata structure contains
%  - colormap = individual colormap, if empty, look for pdata(1), then default
%  - cilm     = individual c-axis, if empty, look for pdata(1), then default [0 N-1].
function FIG = pdata_pcolor(pdata)
FIG.ha = zeros(1, numel(pdata));
FIG.hs = zeros(1, numel(pdata));
FIG.hc = zeros(1, numel(pdata));
FIG.ht = zeros(1, numel(pdata));
xx = pdata(1).xx;
yy = pdata(1).yy;
for k = 1:numel(pdata)
    % Create the plot domain
    FIG.ha(k) = axes('Unit', 'Normalized');
    % Colormap
    if ~isempty(pdata(k).colormap)
        cname = pdata(k).colormap;
        if k == 1
            % Possibly a common colormap case
            colormap(boonlib(cname));
        else
            % Otherwise, we can start individual colormap at k >= 2
            colormap(FIG.ha(k), boonlib(cname));
        end
    elseif ~isempty(pdata(1).colormap)
        % Common colormap case
        cname = pdata(1).colormap;
    end
    % Plot the data
    if exist('cname', 'var') && strcmp(cname(1:3), 'nws') && strcmp(cname(5:7), 'map')
        % The function to convert to index
        ind = boonlib(['nws', cname(4), '2ind'], pdata(k).value);
        FIG.hs(k) = pcolor(xx, yy, ind);
        FIG.hc(k) = boonlib(['nws', cname(4), 'bar']);
    else
        FIG.hs(k) = pcolor(xx, yy, pdata(k).value);
        % CLim
        if ~isempty(pdata(k).clim)
            caxis(pdata(k).clim);
        elseif ~isempty(pdata(1).clim)
            caxis(pdata(1).clim);
        end
        FIG.hc(k) = colorbar;
    end
    % Put the title string
    FIG.ht(k) = title(pdata(k).title);
end
return

function FIG = pdata_arrange(FIG, config)
ndig = ceil(log10(config+1));
if ndig == 3
    nrow = floor(config / 100);
    ncol = rem(floor(config * 0.1), 10);
    tsp = rem(config, 10);
elseif ndig == 2
    nrow = floor(config / 10);
    ncol = rem(config, 10);
    tsp = 0;
else
    error('Unable to understand the parameter %s\n', config);
end

% ox = 0.055;
% oy = 0.056;

set(gcf, 'PaperPositionMode', 'Auto');
tmp = get(gcf, 'Position');
ox = max(0.055, 55 / tmp(3));  % space for y-label
oy = max(0.055, 40 / tmp(4));  % space for x-label

% rect origin
if tsp > 0
    ry = max((tmp(4) - 60) / tmp(4), 0.87) / nrow;
else
    ry = 0.95 / nrow;
end
% rect height
rh = 0.88 * ry;

if all(isempty(FIG.hc(2:end)))
    rx = 0.88 / ncol;  % rect origin
    rw = 0.9 * rx;     % rect width
    for k = 1 : nrow * ncol
        icol = rem(k - 1, ncol);
        irow = floor((k - 1) / ncol);
        set(FIG.ha(k), 'Position', [ox + icol * rx, oy + (nrow - 1 - irow) * ry, rw, rh]);
        if irow < nrow - 1, set(FIG.ha(k), 'XTickLabel', []); else, xlabel(FIG.ha(k), 'X-Distance (m)'); end
        if icol > 0, set(FIG.ha(k), 'YTickLabel', []); else, ylabel(FIG.ha(k), 'Y-Distance (m)'); end
    end
    delete(FIG.hc(2:end))
    set(FIG.hc(1), 'Position', [ox + (ncol - 1) * rx + rw + 0.025, oy + (nrow - 1) * ry, 0.07 * rh, rh])
    FIG.hc = FIG.hc(1);
else
    rx = (1.0 - ox) / ncol;   % rect origin
    rw = 0.80 * rx;           % rect width
    for k = 1 : nrow * ncol
        icol = rem(k - 1, ncol);
        irow = floor((k - 1) / ncol);
        set(FIG.ha(k), 'Position', [ox + icol * rx, oy + (nrow - 1 - irow) * ry, rw, rh]);
        if irow < nrow - 1, set(FIG.ha(k), 'XTickLabel', []); else, xlabel(FIG.ha(k), 'X-Distance (m)'); end
        if icol > 0, set(FIG.ha(k), 'YTickLabel', []); else, ylabel(FIG.ha(k), 'Y-Distance (m)'); end
        set(FIG.hc(k), 'Position', [ox + icol * rx + rw + 0.25 * ox, oy + (nrow - 1 - irow) * ry, 0.03 * rh, rh])
    end
end
if tsp > 0
    FIG.ta = axes('Unit', 'Normalized', 'Position', [0.5 0.94 0.01 0.01]);
    FIG.tt = title('Montage Title', 'FontSize', 14, 'FontWeight', 'Bold');
    axis off
end
return
