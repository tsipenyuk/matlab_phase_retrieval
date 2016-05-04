function plot(obj, varargin)
% PLOT plotting routine for phaseRetrievalAlgorithm objects
%
% Optional arguments: axes handles to axes where the following 
% will be plotted:
%     ax1    solution desity
%     ax2    solution FT modulus and phase
%     ax3    current density
%     ax4    current FT modulus and phase
%     ax5    list of energies of approximations
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
    if obj.current.grid.dimension == 1
        p = inputParser;
        
        addOptional(p,'ax1',0);
        addOptional(p,'ax2',0);
        addOptional(p,'ax3',0);
        addOptional(p,'ax4',0);
        addOptional(p,'ax5',0);
        
        parse(p,varargin{:});
        
        if length(varargin) == 5
            ax1 = p.Results.ax1;
            ax2 = p.Results.ax2;
            ax3 = p.Results.ax3;
            ax4 = p.Results.ax4;
            ax5 = p.Results.ax5;
        elseif length(varargin) == 0
            fig = figure();
            subplot(2,3,1)
            ax1 = gca();
            subplot(2,3,2)
            ax2 = gca();
            subplot(2,3,3)
            ax3 = gca();
            subplot(2,3,4)
            ax4 = gca();
            subplot(2,3,5)
            ax5 = gca();
        else
            disp(strcat(...
                ['You must provide either five optional arguments ' ...
                 '(handles to axes where the ' ...
                 'information will be plotted), ' ...
                 'or none. ' ...
                 'Provided: '], ...
                str(varargin), ', i.e. ', str(length(varargin)), ...
                ' arguments.'))
        end
        
        obj.solution.plot(ax1,ax4);
        title(ax1,'\fontsize{12}Solution density')
        title(ax4,'\fontsize{12}Solution FT')
        obj.current.plot(ax2,ax5);
        title(ax2,'\fontsize{12}Current density')
        title(ax5,'\fontsize{12}Current FT')
        plot(ax3, obj.energy);
        title(ax3,'\fontsize{12}Energy')
    else
        disp(['Currently, only 1-dimensional plotting method ' ...
              'is supported'])
    end
end
