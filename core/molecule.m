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
                continue
            else
                construction_method = 'first_method';
            end

    end
end