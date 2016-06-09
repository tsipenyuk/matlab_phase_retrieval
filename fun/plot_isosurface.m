function plot_isosurface(x1,x2,x3,g)
% plot_isosurface - Custom plot for 3D-densities
%
% Synopsis
%   plot_isosurface(x1,x2,x3,g)
%
% Description
%   Plots transparent isosurfaces corresponding to the density
%   values max(density) * 10.^(-i), i from 1 to 4.
%
% Inputs
%   (meshgrid) x1     First dimension meshgrid
%   (meshgrid) x2     Second dimension meshgrid
%   (meshgrid) x3     Third dimension meshgrid
%   (meshgrid) g      Density values
%
% Examples
%   %% 3D, two gaussians
%   [x1, x2, x3] = meshgrid([-20:0.2:20], [-20:0.2:20],  [-20:0.2:20]);
%   g = exp(-x1.^2 - x2.^2 - x3.^2);
%   shift = fix(length(g)/4);
%   g_sol = circshift(g, [0, shift, shift]) + ...
%           circshift(g, [0, -shift, 0]);
%   plot_isosurface(x1,x2,x3,g_sol);
%
% References
%   http://de.mathworks.com/help/matlab/ref/isosurface.html
%       --- MATLAB documentation on 'isosurface'
%   http://de.mathworks.com/matlabcentral/newsreader/view_thread/268693
%       --- related MATLAB newsreader thread
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase nRetrieval Sandbox root folder.
%
% Changes
%   2016-06-08  First Edition
    figure; 
    hold on; 
    n = 4;
    isolimits = 10.^(1:n);
    colorlist = jet(n);
    colorlist = colorlist(end:-1:1,:);

    for i = 1:1:n
        a = isolimits(i); % 'a' defines the isosurface limits
        p = patch(isosurface(x1,x2,x3,g,max(g(:))/a));
        isonormals(x1,x2,x3,g,p);
        p.FaceColor = colorlist(i,:); % Multicolor isosurfaces
        %p.FaceColor = 'red'; % Monochrome isosurfaces
        p.EdgeColor = 'none';
        daspect([1,1,1]);
        view(3); 
        axis tight;
        alpha(.1)
    end
    
end