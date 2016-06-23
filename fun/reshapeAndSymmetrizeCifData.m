function out = reshapeAndSymmetrizeCifData(varargin)
% reshapeAndSymmetrizeCifData - reshape cifData to a structure of ndgrid arrays
%
% Synopsis 
%   out = reshapeAndSymmetrizeCifData(varargin)
%
% Description
%   Parces a structure returned by the function
%   cif2mat. One may either pass the structure and the necessary
%   field names, or one may pass the value vectors directly. By
%   default, this function assumes that the diffraction data is
%   described in the first octant of the 3d-space, and performs the
%   symmetrization (reflecting the values in the h, k, and l-axes).
%
% Inputs ([]s are optional)
%   Variant 1
%   Necessary arguments:
%   (struct) cifData    Output of the cif2mat function
%   Parameter arguments:
%   Type     Name            Default value  Interpretation
%   (string) 'limBlockName'  'reflns'       Block where max/min
%                                           values of the hkl-grid
%                                           are stored
%   (string) 'hMinName'      'limit_h_min'
%   (string) 'hMaxName'      'limit_h_max' 
%   (string) 'kMinName'      'limit_k_min' 
%   (string) 'kMaxName'      'limit_k_max' 
%   (string) 'lMinName'      'limit_l_min' 
%   (string) 'lMaxName'      'limit_l_max' 
%   (string) 'valBlockName'  'refln'        Block where indices and
%                                           measurement values are stored
%   (string) 'hName'         'index_h' 
%   (string) 'kName'         'index_k' 
%   (string) 'lName'         'index_l' 
%   (string) 'fName'         'F_meas_au'    Intensity values
%   (string) 'fSigmaName'    'None'         Intensity uncertainty
%   (bool)   'sym'           true           Should the data be symmetrized?
%
%   Variant 2
%   Parameter arguments:
%   Type     Name            Default value  Interpretation
%   (float)  'hMin'                         Min of h-axis
%   (float)  'hMax'                         Max of h-axis
%   (float)  'kMin'                         Min of k-axis
%   (float)  'kMax'                         Max of k-axis
%   (float)  'lMin'                         Min of l-axis
%   (float)  'lMax'                         Max of l-axis
%   (vector) 'h'                            Values of h
%   (vector) 'k'                            Values of k
%   (vector) 'l'                            Values of l
%   (vector) 'f'                            Values of intensity
%   (vector) 'fSigma'        'None'         Values of uncertainty (optional)
%   (bool)   'sym'           true           Should the data be
%                                           symmetrized? (optional)
%
% Outputs ([]'s are optional)
%   (struct) out  structure containing the following fields:
%   Field name    Contents
%   out.h         Same as input, vector containing h-entries
%   out.k         Same as input, vector containing k-entries
%   out.l         Same as input, vector containing l-entries
%   out.f         Same as input, vector containing f-entries
%   [out.fSigma]  Same as input, vector containing uncertainty entries
%   out.H         ndarray, box from hMin to hMax, from kMin to
%                 kMax, from lMin to lMax;
%   out.K         ndarray, --//--
%   out.L         ndarray, --//--
%   out.F         ndarray with values f(i) at h(i), k(i), l(i),
%                 padded with zeros where F is not specified
%   [out.FSigma]  --//-- uncertainties fSigma(i)
%
%
% Examples
%  % Example 1
%  >> cifData = cif2mat('data/2OLX/2olx-sf.cif')
%  >> rData = reshapeAndSymmetrizeCifData(cifData);
%  >> sqrtI = rData.F;
%  >> figure; h = slice(permute(rData.H,[2,1,3]),...
%  permute(rData.K,[2,1,3]), permute(rData.L,[2,1,3]),...
%  permute(sqrtI,[2,1,3]), 0, 0, 0);
%   
%  % Example 2
%  >> cifData = cif2mat('data/5B3D/5b3d-sf.cif')
%  >> rData = reshapeAndSymmetrizeCifData(cifData, 'limBlockName',...
%  'diffrn_reflns', 'fName', 'intensity_meas');
%  >> sqrtI = rData.F;
%  >> figure; h = slice(permute(rData.H,[2,1,3]),...
%  permute(rData.K,[2,1,3]), permute(rData.L,[2,1,3]),...
%  permute(sqrtI,[2,1,3]), 0, 0, 0);
%
% See also
%   cif2mat
%   
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase nRetrieval Sandbox root folder.
%
% Changes
%   2016-06-23  First Edition
     
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
        p.addParamValue('limBlockName', 'reflns', @isstr);
        p.addParamValue('hMinName', 'limit_h_min', @isstr);
        p.addParamValue('hMaxName', 'limit_h_max', @isstr);
        p.addParamValue('kMinName', 'limit_k_min', @isstr);
        p.addParamValue('kMaxName', 'limit_k_max', @isstr);
        p.addParamValue('lMinName', 'limit_l_min', @isstr);
        p.addParamValue('lMaxName', 'limit_l_max', @isstr);
        p.addParamValue('valBlockName', 'refln', @isstr);
        p.addParamValue('hName', 'index_h', @isstr);
        p.addParamValue('kName', 'index_k', @isstr);
        p.addParamValue('lName', 'index_l', @isstr);
        p.addParamValue('fName', 'F_meas_au', @isstr);
        p.addParamValue('fSigmaName', 'None', @isstr);
        p.addParamValue('sym', true, @islogical);
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
        p.addParamValue('fSigma', 'None', @isvector);
        p.addParamValue('sym', true, @islogical);
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
                out.FSigma( out.h(iVal) + h0,  ...
                            out.k(iVal) + k0,  out.l(iVal) + l0) ...
                    = out.fSigma(iVal);
            end
    end
    
end