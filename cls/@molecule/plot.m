function plot(obj, varargin)
% PLOT plotting routine for molecule density and its FT
%
% Optional arguments: axes handles to axes where the
% density and the FT modulus and phase will be plotted
%
% Examples
% Simple plot
%   >> m = molecule()
%   >> m.plot()
% 
% Plot using specific axes
%   >> figure();
%   >> subplot(2,2,1)
%   >> ax1 = gca
%   >> subplot(2,2,4)
%   >> ax2 = gca
%   >> m.plot(ax1, ax2)
    if obj.grid.dimension == 1
        p = inputParser;
        
        addOptional(p,'ax1',0);
        addOptional(p,'ax2',0);
        
        parse(p,varargin{:});
        
        if length(varargin) == 2
            ax1 = p.Results.ax1;
            ax2 = p.Results.ax2;
        elseif length(varargin) == 0
            fig = figure();
            subplot(2,1,1)
            ax1 = gca();
            subplot(2,1,2)
            ax2 = gca();
        else
            disp(strcat(...
                ['You must provide either two optional arguments ' ...
                 '(ax1, ax2, handles to axes where the ' ...
                 'density and FT  will be plotted), ' ...
                 'or none. ' ...
                 'Provided: '], ...
                str(varargin), ', i.e. ', str(length(varargin)), ...
                ' arguments.'))
        end
        
        hold on
        % 1st axis
        plot(ax1, obj.grid.xAxes{1}, obj.density)
        try % If density == const, this may cause problems
            axis(ax1, [min(obj.grid.xAxes{1}) ...
                       max(obj.grid.xAxes{1}) ...
                       min(obj.density)...
                       max(obj.density)])
            %ax1.XTick = [0 pi/2 pi  3*pi/2 2*pi];
            %ax1.XTickLabel = {'0', '\pi/2', '\pi',  '3\pi/2', ...
            %                  '2\pi'};
        end
        title(ax1,'\fontsize{12}Density')
        
        
        
        % 2nd axis
        yyaxis(ax2, 'left')
        plot(ax2, obj.grid.kAxes{1}, abs(fftshift(obj.ft)))
        try
            axis(ax2, [min(obj.grid.kAxes{1}) ...
                       max(obj.grid.kAxes{1}) ...
                       min(abs(fftshift(obj.ft))) ...
                       max(abs(fftshift(obj.ft)))])
        end
        title(ax2,'\fontsize{12}FT: modulus and phase')
        
        % 3nd axis
        yyaxis(ax2, 'right')
        plot(ax2, obj.grid.kAxes{1}, ...
             angle(fftshift(obj.ft)))
        try
            axis(ax2, [min(obj.grid.kAxes{1}) ...
                       max(obj.grid.kAxes{1}) ...
                       -pi ...
                       pi])
        end
        hold off
        
    elseif obj.grid.dimension == 2
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
        image(ax1, obj.grid.xAxes{1}, obj.grid.xAxes{2}, ...
              obj.density, 'CDataMapping', 'scaled')
        %axis equal tight
        title(ax1,'\fontsize{12}Density')
        colormap(linear_bgyw_15_100_c68_n256)
        colorbar(ax1)
        
        % 1st axis
        image(ax2, obj.grid.kAxes{1}, obj.grid.kAxes{2}, ...
              abs(fftshift(obj.ft)), 'CDataMapping', 'scaled')
        %axis equal tight
        title(ax2,'\fontsize{12}FT modulus')
        colormap(linear_bgyw_15_100_c68_n256)
        colorbar(ax2)
        
        % 3nd axis
        image(ax3, obj.grid.kAxes{1}, obj.grid.kAxes{2}, ...
              angle(fftshift(obj.ft)), 'CDataMapping', 'scaled')
        %axis equal tight                
        title(ax3,'\fontsize{12}FT phase')
        colormap(ax3, cyclic_mygbm_30_95_c78_n256)
        colorbar(ax3)
        
    else
    disp(['Currently, only 1-dimensional plotting method ' ...
          'is supported'])
    end
end
