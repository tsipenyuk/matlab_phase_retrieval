function contour(obj, varargin)
% CONTOUR plotting routine for 2-dim molecule density and its FT
%
% Optional arguments: axes handles to axes where the
% density and the FT modulus and phase will be plotted
    if obj.grid.dimension == 2
        p = inputParser;

        addOptional(p,'ax1',0);
        addOptional(p,'ax2',0);
        addOptional(p,'ax3',0);

        parse(p,varargin{:});

        if length(varargin) == 3
            ax1 = p.Results.ax1;
            ax2 = p.Results.ax2;
            ax3 = p.Results.ax3;
        elseif length(varargin) == 0
            fig = figure();
            subplot(2,2,1)
            ax1 = gca();
            subplot(2,2,2)
            ax2 = gca();
            subplot(2,2,3)
            ax3 = gca();
        else
            disp(strcat(...
                ['You must provide either three optional arguments ' ...
                 '(ax1, ax2, ax3 handles to axes where the ' ...
                 'density and FT modulus and phase' ...
                 '  will be plotted), or none. ' ...
                 'Provided: '], ...
                str(varargin), ', i.e. ', str(length(varargin)), ...
                ' arguments.'))
        end


        colormap linear_bgyw_15_100_c68_n256

        % 1st axis
        contour(ax1, obj.grid.xAxes{1}, obj.grid.xAxes{2}, ...
                obj.density)
        %axis equal tight
        title(ax1,'\fontsize{12}Density')
        colormap(linear_bgyw_15_100_c68_n256)
        colorbar(ax1)

        % 1st axis
        contour(ax2, obj.grid.kAxes{1}, obj.grid.kAxes{2}, ...
                abs(fftshift(obj.ft)))
        %axis equal tight
        title(ax2,'\fontsize{12}FT modulus')
        colormap(linear_bgyw_15_100_c68_n256)
        colorbar(ax2)

        % 3nd axis
        contour(ax3, obj.grid.kAxes{1}, obj.grid.kAxes{2}, ...
                angle(fftshift(obj.ft)))
        %axis equal tight                
        title(ax3,'\fontsize{12}FT phase')
        colormap(ax3, cyclic_mygbm_30_95_c78_n256)
        colorbar(ax3)

    else
        disp(['Currently, only 2-dimensional contour plot ' ...
              'is supported'])
    end
end
