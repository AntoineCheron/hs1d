function [x,y] = ind2coord(FD,ix)

% convert linear index to x and y coordinates
%
% Syntax
%
%     [x,y] = ind2coord(FD,ix)
%
% Description
%
%     ind2coord converts a linear index into an instance of FLOWobj to x
%     and y coordinates.
%
% Input arguments
%
%     FD      instance of FLOWobj
%     ix      linear index
%
% Output arguments
%
%     x,y     x- and y-coordinates  
%
% See also: FLOWobj/coord2ind
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 08. June, 2015


[r,c] = ind2sub(FD.size,ix(:));
xy    = double([r c ones(numel(ix),1)])*FD.refmat;
x = xy(:,1);
y = xy(:,2);
