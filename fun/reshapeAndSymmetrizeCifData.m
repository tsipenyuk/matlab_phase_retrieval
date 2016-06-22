function out = reshapeAndSymmetrizeCifData(varargin)
% reshapeCifData - reshape cifData to a structure of ndgrid arrays
%
% Synopsis 
%   [g_new, error] = dmap_optimized(g, sqrtI, [beta], [pObj], [varargin])
%
% Description
%   P
%
% Inputs ([]s are optional)
%   (ndarray) g      Current  approximation to the solution of the
%                    phase retrieval problem
%   (ndarray) sqrtI      Phase retrieval data (square root of the
%                    measured intensity)
%   (scalar)  [beta = 1]
%                    Update parameter, default value reproduces the
%                    hybrid input-output update (Bauschke' variant)
%   (func)    [pObj = @pP]
%                    handle to the projection onto the object space
%                    (non-negativity, support, etc.). pObj must
%                    take g as the first argument, may contain
%                    optional arguments specified in varargin.
%   (...)     [varargin]
%                    optional arguments passed to pObj --- if
%                    submitted, pObj is  called as pObj(g, varargin)
%
% Outputs ([]s are optional)
%   (ndarray) g_new  Updated approximation to the solution of the 
%                    phase retrieval problem
%   (scalar)  error  Error (energy) corresponding to the updated 
%                    approximation
%
% Examples
%   %% 1D, two gaussians
%   x = [-20:0.2:20];
%   g_sol = exp(-x.^2);
%   shift = fix(length(g)/4);
%   g_sol = circshift(g, [0, shift]) + circshift(g, [0, -shift]);
%   sqrtI = abs(fftn(g_sol));
%   g_new = zeros(size(g));
%   E = [];
%   for i=1:1:200
%       [g_new, error] = dmap_optimized(g_new, sqrtI);
%       E = [E error];
%   end
%   plot(E);
%   plot(g_new);
%
%   %% 2D, two gaussians
%   [x1, x2] = ndgrid([-20:0.2:20], [-20:0.2:20]);
%   g = exp(-x1.^2 - x2.^2);
%   shift = fix(length(g)/4);
%   g_sol = circshift(g, [0, shift]) + circshift(g, [0, -shift]);
%   sqrtI = abs(fftn(g_sol));
%   g_new = zeros(size(g));
%   E = [];
%   for i=1:1:200
%       [g_new, error] = dmap_optimized(g_new, sqrtI);
%       E = [E error];
%   end
%
% See also
%   er
%   bio
%   hio_fienup
%   dmap
%   
% Requirements
%   pM (modulus projection)
%   pP (non-negative projection)
%
% References
%   V. Elser, “Phase retrieval by iterated projections,” 
%   J. Opt. Soc. Am. A., vol. 20, pp. 40–55, 2003.
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase nRetrieval Sandbox root folder.
%
% Changes
%   2016-06-10  First Edition
     
    % Check which data type was passed
    iscifDataPassed = false;
    try
        if varargin{1}.type == 'cifData';
            iscifDataPassed = true;
        end
    end
    
    % Parse input
    p = inputParser;
    if iscifDataPassed == true
        p.addRequired('cifData', @isstruct);
        p.addOptional('limBlockName', 'diffrn_reflns', @isstr);
        p.addOptional('hMinName', 'limit_h_min', @isstr);
        p.addOptional('hMaxName', 'limit_h_max', @isstr);
        p.addOptional('kMinName', 'limit_k_min', @isstr);
        p.addOptional('kMaxName', 'limit_k_max', @isstr);
        p.addOptional('lMinName', 'limit_l_min', @isstr);
        p.addOptional('lMaxName', 'limit_l_max', @isstr);
        p.addOptional('valBlockName', 'refln', @isstr);
        p.addOptional('hName', 'index_h', @isstr);
        p.addOptional('kName', 'index_k', @isstr);
        p.addOptional('lName', 'index_l', @isstr);
        p.addOptional('fName', 'intensity_meas', @isstr);
        p.addOptional('fSigmaName', 'None', @isstr);
        p.addOptional('sym', true, @islogical);
    else
        p.addRequired('hMin', @isfloat);
        p.addRequired('hMax', @isfloat);
        p.addRequired('kMin', @isfloat);
        p.addRequired('kMax', @isfloat);
        p.addRequired('lMin', @isfloat);
        p.addRequired('lMax', @isfloat);
        p.addRequired('h', @isvector);
        p.addRequired('k', @isvector);
        p.addRequired('l', @isvector);
        p.addRequired('f', @isvector);
        p.addOptional('fSigma', 'None', @isvector);
        p.addOptional('sym', true, @islogical);
    end
    
    p.parse(varargin{:});
    pR = p.Results;
    if iscifDataPassed == true
        hMin = pR.cifData.(pR.limBlockName).(pR.hMinName);
        hMax = pR.cifData.(pR.limBlockName).(pR.hMaxName);
        kMin = pR.cifData.(pR.limBlockName).(pR.kMinName);
        kMax = pR.cifData.(pR.limBlockName).(pR.kMaxName);
        lMin = pR.cifData.(pR.limBlockName).(pR.lMinName);
        lMax = pR.cifData.(pR.limBlockName).(pR.lMaxName);
        h = pR.cifData.(pR.valBlockName).(pR.hName);
        k = pR.cifData.(pR.valBlockName).(pR.kName);
        l = pR.cifData.(pR.valBlockName).(pR.lName);
        f = pR.cifData.(pR.valBlockName).(pR.fName);
        if strcmp(pR.fSigmaName, 'None') == true
            fSigma = 'None';
        else
            fSigma = pR.(valBlockName).(fSigmaName);
        end
    else
        hMin = pR.hMin;
        hMax = pR.hMax;
        kMin = pR.kMin;
        kMax = pR.kMax;
        lMin = pR.lMin;
        lMax = pR.lMax;
        h = pR.h;
        k = pR.k;
        l = pR.l;
        f = pR.f;
        fSigma = pR.fSigma;
    end
    
    % keep the old info - to be able to access it with unified notation
    out.h = h;
    out.k = k;
    out.l = l;
    out.f = f;
    out.fSigma = fSigma;
    h0 = 1; % Index of the coordinates origin
    k0 = 1; 
    l0 = 1; 

    % If needed, the data is symmetrized, assuming that the
    % provided data lies in the top-front-right (+ + +) octant
    if pR.sym == true
        disp('dang')
        out.h = [ h -h -h  h  h -h -h  h];
        out.k = [ k  k -k -k  k  k -k -k];
        out.l = [ l  l  l  l -l -l -l -l];
        out.f = [ f  f  f  f  f  f  f  f];
        if strcmp(fSigma, 'None') ~= true
            out.fSigma = [ fSigma  fSigma  fSigma  fSigma ...
                           fSigma  fSigma  fSigma  fSigma];
        end
        
        hMin = -hMax;
        kMin = -kMax;
        lMin = -lMax;
        h0 = hMax + 1; % Index of the coordinates origin
        k0 = kMax + 1; 
        l0 = lMax + 1; 
    end

    [out.H, out.K, out.L] = ndgrid(hMin:hMax, kMin:kMax, lMin:lMax);
    out.F = zeros(size(out.H));

    disp(length(out.h))
    disp(length(out.k))
    disp(length(out.l))
    disp(length(out.f))
    for iVal = 1:1:length(out.f)
        out.F(out.h(iVal) + h0, out.k(iVal) + k0, out.l(iVal) + l0) = out.f(iVal);
    end
    
    if strcmp(fSigma, 'None') ~= true
            for iVal = 1:1:length(out.f)
                out.F_sigma( out.h(iVal) + h0,  ...
                             out.k(iVal) + k0,  out.l(iVal) + l0) ...
                    = out.f_sigma(iVal);
            end
    end
    
end