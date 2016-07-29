function plot_fnd(zGrid, fGrid, gGrid, fName, gName, ax, varargin)
% plot_fnd - Quiver plots of two complex-valued functions
%
% Synopsis ([]'s are optional)
%   plot_isosurface(xGrid, fHandle, dfHandle, [fName], [gName], varargin)
%
% Description
%   Quiver plot of functions fGrid and gGrid from $|mathbb C$
%   to  $\mathbb C$ evaluated at zGrid, 2-dimensional array of
%   complex points.
%
% Inputs
%   (grid)   xGrid          2-dimensional grid of complex numbers
%   (grid)   fGrid          2-dimensional grid of function f values
%   (grid)   gGrid          2-dimensional grid of function g values
%   (string) [fName='fun1'] function label in the plot
%   (string) [gName='fun2'] function label in the plot
%
% Examples
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase nRetrieval Sandbox root folder.
%
% Changes
%   2016-07-26  First Edition

if nargin == 3
    fName = 'fun1';
    gName = 'fun2';
end

if nargin < 6
    figure;
    ax = axes();
end

hold(ax, 'on');
%quiver(ax, real(zGrid), imag(zGrid), real(fGrid), imag(fGrid), ...
%       varargin{:});
contour(ax,real(zGrid), imag(zGrid), real(fGrid))
quiver(ax, real(zGrid), imag(zGrid), real(gGrid), imag(gGrid), varargin{:});
legend(ax, fName, gName,'Interpreter','LaTex')
axis tight
axis equal
hold(ax, 'off');
end
