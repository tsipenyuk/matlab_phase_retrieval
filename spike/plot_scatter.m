function plot_scatter(h,k,l,f)
% plot_scatter - 3d scatter plot for intensity measurements
%
% Synopsis
%   plot_scatter(h,k,l,f)
%
% Description
%   Scatter plot, plots points at coordinates h,k,l with intensity
%   f 
%
% Inputs
%   (vector) h     First dimension coordinates
%   (vector) k     Second dimension coordinates
%   (vector) l     Third dimension coordinates
%   (vector) f     Scattering intensity
%
% Examples
%   %% 3D gaussian
%   [H, K, L] = ndgrid([-2:0.5:2], [-2:0.5:2], [-2:0.5:2]);
%   h = reshape(H, [], 1);
%   k = reshape(K, [], 1);
%   l = reshape(L, [], 1);
%   f = exp(-h.^2 - k.^2 - l.^2);
%   plot_scatter(h,k,l,f);
%
% References
%   http://de.mathworks.com/help/matlab/ref/scatter3.html
%       --- MATLAB documentation on 'scatter3'
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase nRetrieval Sandbox root folder.
%
% Changes
%   2016-06-15  First Edition
    figure; 
    colormap(linear_kryw_5_100_c67_n256); 
    colormap(flipud(colormap));
    edgeColor = [0.8 0.8 0.8];
    scatter3(h, k, l, 30, f, 'filled', ...
             'MarkerEdgeColor', edgeColor); 
    colorbar;
end