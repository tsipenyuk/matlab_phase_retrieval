function plot_slice(H,K,L,F,Sh,Sk,Sl,gridtype, varargin)
% plot_slice - Custom slice plot for intensity measurements
%
% Synopsis ([]'s are optional)
%   plot_isosurface(H,K,L,F,Sh,Sk,Sl,[gridtype])
%
% Description
%   Slice plot for intensity mesurements. The first four arguments
%   must be all either ndgrids (by default) or meshgrids. In the
%   latter case, set 'gridtype' to 'm'.
%
% Inputs
%   (grid) H     First dimension meshgrid or ndgrid
%   (grid) K     Second dimension meshgrid or ndgrid
%   (grid) L     Third dimension meshgrid or ndgrid
%   (grid) F     Intensity values, meshgrid or ndgrid
%
% Examples
%   %% 3D, two gaussians
%   [H, K, L] = ndgrid([-20:0.2:20], [-20:0.2:20],  [-20:0.2:20]);
%   g = exp(-x1.^2 - x2.^2 - x3.^2);
%   shift = fix(length(g)/4);
%   F = circshift(g, [0, shift, shift]) + ...
%       circshift(g, [0, -shift, 0]);
%   plot_slice(H,K,L,F,[-1.2 .8 2],2,[-2 -.2])
%
% References
%   http://de.mathworks.com/help/matlab/ref/slice.html
%       --- MATLAB documentation on 'slice'
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
    if nargin < 8
        gridtype = 'n';
    end
    
    if nargin < 7
        F = H;
        Sh = [1 20];
        Sk = 1;
        Sl = 1;
        if gridtype == 'n'
            [H, K, L] = ndgrid(1:size(F,1),1:size(F,2),1:size(F,3));
        else
            [H, K, L] = meshgrid(1:size(F,1),1:size(F,2),1:size(F, ...
                                                              3));
        end
    end
    disp(size(H))
    disp(size(K))
    disp(size(L))
    
    colormap(linear_kryw_5_100_c67_n256); 
    colormap(flipud(colormap));

    if gridtype == 'm'
        h = slice(H, K, L, F, Sh, Sk, Sl, varargin);
    else
        h = slice(permute(H,[2,1,3]), permute(K,[2,1,3]), ...
                  permute(L,[2,1,3]), permute(F,[2,1,3]), Sh, ...
                  Sk, Sl, varargin);
    end
    edgeColor = [0.8 0.8 0.8];
    set(h, 'EdgeColor', edgeColor);
    colorbar;
end