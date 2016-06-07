function f_list = rptCoordinates(x_list, x_com)
% RPTCOORDINATES Given center-of-mass coordinates of mass pieces, calculate
% their right-point coordinates assuming constant weight
% distribution within one piece.
%
%   x_rpt = RPTCOORDINATES(x_pcs, f_pcs) Given...
%
%   Arseniy Tsipenyuk, TUM M7
%   May 10th, 2016
    
    x_list = [left_bd zeros(size(x_com)) right_bd];
    x_rpt = (x_list(2:end) + x_list(1:end-1)) / 2;
end
