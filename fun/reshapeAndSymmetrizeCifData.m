function out = reshapeAndSymmetrizeCifData(varargin)
% reshapeAndSymmetrizeCifData - reshape cifData to a structure of ndgrid arrays
%
% Synopsis 
%   out = reshapeAndSymmetrizeCifData(varargin)
%
% Examples
%  % Example 1
%  cifData = cif2mat('data/2y29/2y29-sf.cif')
%  rData = reshapeAndSymmetrizeCifData(cifData, 'fName', 'F_meas_au', 'padZeros',2);
%  figure; h = slice(permute(rData.H,[2,1,3]),...
%     permute(rData.K,[2,1,3]), permute(rData.L,[2,1,3]),...
%     permute(sqrtI,[2,1,3]), 0, 0, 0);
%  sqrtI = rData.F;
%  sqrtIpadded = ones(size(sqrtI));
%  sqrtIpadded(~isnan(sqrtI)) = sqrtI(~isnan(sqrtI));
%  gInit = pP(real(ifftn(sqrtIpadded.*sqrtIpadded)));
%  [gOut eOut] = compare_algorithms(gInit, sqrtI, 10000, @er, ...
%  @hio_fienup);
%  figure; hist(gOut{1}(:), 10)
%  [X, Y, Z] = meshgrid(1:size(gOut{1},1),1:size(gOut{1},2), 1:size(gOut{1},3));
%  % X = X * rData.len_a; Y = Y * rData.len_b; Z = Z * rData.len_c;
%  gS = smooth3(gOut{1});
%  figure; hist(gS(:),10);
%  figure; isosurface(X, Y, Z, permute(smooth3(gOut{1}),[2,1,3]),
%  0.32); axis equal;
%
% Description
%   Parces a structure returned by the function
%   cif2mat. One may either pass the structure and the necessary
%   field names, or one may pass the value vectors directly. By
%   default, this function assumes that the diffraction data is
%   described in the first octant of the 3d-space, and performs the
%   symmetrization.
%
% Caveat
%   Does not respect slanted or transformed cell structure. In
%   other words, if the angles between cell axes are not equal to
%   90 degrees or if a cell transformation took place
%   (cf. _diffrn_reflns_transf_matrix_), the returned coordinates
%   of out.h, out.k, out.l, out.H, out.K, out.L are not correct. 
%
% Inputs ([]s are optional)
%   Variant 1
%   Necessary arguments:
%   (struct) cifData    Output of the cif2mat function
%   Parameter arguments:
%   Type     Name            Default value  Interpretation
%   (string) 'cellBlockName' 'cell'         Block where cell lengths
%                                           are stored
%   (string) 'aName'         'length_a'
%   (string) 'bName'         'length_b'
%   (string) 'cName'         'length_c'
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
%   out.h         Vector containing h-entries, scaled by cell length
%   out.k         Vector containing k-entries, scaled by cell length
%   out.l         Vector containing l-entries, scaled by cell length
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
        p.addParamValue('cellBlockName', 'cell', @isstr);
        p.addParamValue('aName', 'length_a', @isstr);
        p.addParamValue('bName', 'length_b', @isstr);
        p.addParamValue('cName', 'length_c', @isstr);
        p.addParamValue('valBlockName', 'refln', @isstr);
        p.addParamValue('hName', 'index_h', @isstr);
        p.addParamValue('kName', 'index_k', @isstr);
        p.addParamValue('lName', 'index_l', @isstr);
        p.addParamValue('fName', 'F_meas_au', @isstr);
        p.addParamValue('fSigmaName', 'None', @isstr);
        p.addParamValue('sym', true, @islogical);
        p.addParamValue('padZeros', 1, @isfloat);
    else
        p.addRequired('h', @isvector);
        p.addRequired('k', @isvector);
        p.addRequired('l', @isvector);
        p.addRequired('f', @isvector);
        p.addParamValue('fSigma', 'None', @isvector);
        p.addParamValue('sym', true, @islogical);
        p.addParamValue('padZeros', 1, @isfloat);
    end
    
    p.parse(varargin{:});
    pR = p.Results;
    if iscifDataPassed == true
        len_a = pR.cifData.(pR.cellBlockName).(pR.aName);
        len_b = pR.cifData.(pR.cellBlockName).(pR.bName);
        len_c = pR.cifData.(pR.cellBlockName).(pR.cName);
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
        h = pR.h;
        k = pR.k;
        l = pR.l;
        f = pR.f;
        fSigma = pR.fSigma;
        len_a = 1;
        len_b = 1;
        len_c = 1;
    end

    
    hMin = min(h);
    hMax = max(h);
    kMin = min(k);
    kMax = max(k);
    lMin = min(l);
    lMax = max(l);
    % Entries outside this radius will be set to zero
    max_rad = max(max(max(sqrt((h/len_a) .^ 2 +  (k/len_b) .^2 + (l/len_c) .^2))));
    
    
    % Add additional zeros outside of the box 
    % (Oversample to enforce better resolution)
    new_pt_h = hMax * pR.padZeros + 1; % Add 1 to ensure that none
                                       % of the real measurements
                                       % are overridden
    new_pt_k = kMax * pR.padZeros;
    new_pt_l = lMax * pR.padZeros;
    new_pt_f = 0;
    h = [h new_pt_h];
    k = [k new_pt_k];
    l = [l new_pt_l];
    f = [f new_pt_f];
    addedPoints = 1; % One 'fake' point added te ensure oversampling

    % Correct values after oversampling
    hMin = min(h);
    hMax = max(h);
    kMin = min(k);
    kMax = max(k);
    lMin = min(l);
    lMax = max(l);

    
    % keep the old info - to be able to access it with unified notation
    out.h = h;
    out.k = k;
    out.l = l;
    out.f = f;
    out.fSigma = fSigma;
    out.len_a = len_a;
    out.len_b = len_b;
    out.len_c = len_c;
    h0 = 1; % Index of the coordinates origin
    k0 = 1; 
    l0 = 1; 

    % If needed, the data is symmetrized, assuming that the
    % provided data lies in the top-front-right (+ + +) octant
    if pR.sym == true
        out.h = [ h -h -h  h  h -h -h  h];
        out.k = [ k  k -k -k  k  k -k -k];
        out.l = [ l  l  l  l -l -l -l -l];
        out.f = [ f  f  f  f  f  f  f  f];
        addedPoints = 8; % 8 fake points, cf (**) above, were added.
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
    out.F = NaN * zeros(size(out.H));

    
    disp(['Number of diffraction points (after symmetrization): ' num2str(length(out.h)-8) '.'])
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
    
    % Rescale coordinates
    out.h = out.h * 1 / len_a;
    out.k = out.k * 1 / len_b;
    out.l = out.l * 1 / len_c;
    out.H = out.H * 1 / len_a;
    out.K = out.K * 1 / len_b;
    out.L = out.L * 1 / len_c;

    
    % set outer intenisty NaN-values to zero
    zero_val_pts = (sqrt(out.H .^ 2 +  out.K .^2 + out.L .^2) >= ...
                    max_rad);
    out.F(zero_val_pts) = 0;




end