function h = ribboncoloredZ(varargin)
% ribboncoloredZ is a wrapper of buitlin ribbon() and color ribbons
% according to Z values and axes' colormap.
%
% h = ribboncoloredZ(y)
% h = ribboncoloredZ(y,z)
% h = ribboncoloredZ(y,z,width)
% h = ribboncoloredZ(axh,_____)
%
% INPUT ARGUMENTS
% y,z         y can be a matrix whose columns are plotted as separate
%             ribbons. Or, when y and z are provided, they are vectors of
%             the same size.
%
% width       (default) 0.75
%             (Optional) Width of ribbons in X axis
%
% OUTPUT ARGUMENTS
% h           Graphic handle vector of surface objects
%
% See also
% ribbon, surface, colormap, ribboncoloredZ_script
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 29-Jun-2016 03:57:30


narginchk(1,4);
if isscalar(varargin{1}) && isgraphics(varargin{1}) && strcmpi(varargin{1}.Type,'axes')
    assert(nargin >= 2)
    axh = varargin{1};
    axes(axh);
    varargin = varargin(2:end);
else
    figure;
    axh = axes;
end

switch length(varargin)
    case 1
        y = varargin{1};
        h = ribbon(axh,y);
        
    case 2
        y = varargin{1};
        z = varargin{2};
        h = ribbon(axh,y,z);
        
    case 3
        y = varargin{1};
        z = varargin{2};
        width = varargin{3};
        
        h = ribbon(axh,y,z,width);
    otherwise
        error('too many input arguments')
end


for i = 1:length(h)
    h(i).CData = h(i).ZData;
    h(i).FaceColor = 'interp';
    h(i).FaceLighting = 'gouraud';
    h(i).MeshStyle = 'column';
end






end