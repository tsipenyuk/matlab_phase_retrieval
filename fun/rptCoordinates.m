function x_rpt = rptCoordinates(left_bd, x_com)
% RPTCOORDINATES Given center-of-mass coordinates of mass pieces, calculate
% their center-of-mass coordinates assuming linear weight
% distribution within one piece.
    x_rpt_long = horzcat(left_bd, zeros(1, length(x_com)));
    for i = 1 : 1 : length(x_com)
        x_rpt_long(i+1) = 2 * x_com(i) - x_rpt_long(i);
    end
    x_rpt = x_rpt_long(2:end)
end
