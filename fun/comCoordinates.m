function x_com = comCoordinates(left_boundary, x_pcs)
% COMCOORDINATES Given rightmost endpoints of mass pieces, calculate
% their center-of-mass coordinates assuming linear weight
% distribution within one piece.
    x_com = zeros(1, length(x_pcs));
    mass_bounds = horzcat(left_boundary, x_pcs);
    for i = 1 : 1 : length(x_pcs)
        x_com(i) = (mass_bounds(i + 1) + mass_bounds(i)) / 2;
    end
end
