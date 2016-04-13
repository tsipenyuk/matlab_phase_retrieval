classdef phaseRetrievalAlgorithm
%  phaseRetrievalAlgorithm - class to store appoximations and error data
%
% Class attributes:
%   obj.initial  - instance of molecule, stores initial
%                  approximation
%   obj.current  - instance of molecule, stores current
%                  approximation
%   obj.solution - instance of molecule, stores ideal result
%                  (if known)
%   obj.energy   - array storing energies corresponding to all
%                  approximations
%   obj.num_steps -total number of already performed steps
%
% Class methods:
%   obj.update   - perform an update step in the algorithm. 
%
% For examples, type 'doc phaseRetrievalAlgorithm' and see
% Constructor info.
% <UNFINISHED>
    properties
        initial % Initial approximation
        current % Current approximation
    end
    properties (SetAccess=private)
        energy    % Energy vector
        num_steps % Total number of performed steps
        update    % Handle to the update function
    end
    properties (SetAccess=private)
        solution  % Exact solution (if known)
    end
    methods
        function obj = phaseRetrievalAlgorithm();
    end
end