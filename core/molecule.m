classdef molecule
% MOLECULE describes a molecule (its density and its Fourier transform).
%
% <NOT FINISHED>
    properties
        grid
        density
        ft
    end
    methods
        function obj=molecule(varargin)
        % MOLECULE Creates a molecule from density or fourier transform.
        %
        % 
            p = inputParser;

            default_grid = 0; % Dummy default
            default_density = exp(-(linspace(1, 10, 50) - ...
                                    linspace(5, 5, 50)).^2)
            default_ft = fftn(default_density)
            
            addParameter(p,'grid',default_nPoints)
            addParameter(p,'density',default_xMin,@isnumeric)
            addParameter(p,'ft',default_xMax,@isnumeric)

            parse(p,varargin{:});

            % Figure out which arguments were passed and which
            % construction method will be used
            str_varargin = varargin(1:2:size(varargin,2));
            if any(strmatch('density',str_varargin))
                obj.density = p.Results.density
                obj.ft = fftn(obj.density)
            elseif any(strmatch('ft',str_varargin))
                obj.ft = p.Results.ft
                obj.density = real(ifftn(obj.ft))
            else
                disp(['Input density or FT were not specified. ' ...
                      'Proceeding with default density (Gaussian, ' ...
                      'width 1, mean 5, plotted between 0 and ' ...
                      '10).'])
                obj.density = p.Results.density
                obj.ft = p.Results.ft
            end

            if any(strmatch('grid',str_varargin))
                obj.grid = p.Results.grid
            else
                disp('(Warning: grid was not specified.)')
            end
        end

        
        function set_density(obj, density)
            obj.density = density
            obj.ft = fftn(obj.density)
        end

        
        function set_ft(obj, ft)
            obj.ft = ft
            obj.density = real(ifftn(obj.ft))
        end

        
        function set_phase(obj, phase)
            obj.ft = abs(obj.ft) .* exp(1j * phase)
            obj.density = real(ifftn(obj.ft))
        end
        
        
        function set_random_phase(obj)
            random_phase = 2 * pi * rand(size(obj.ft))
            obj.set_phase(random_phase)
        end
        
        
        function set_zero_phase(obj, phase)
            obj.ft = abs(obj.ft)
            obj.density = real(ifftn(obj.ft))
        end

    end
end