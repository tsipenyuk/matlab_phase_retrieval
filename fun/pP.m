function res = pP(g)
% pP Replace negative entries of g by 0
    res = real(abs(g .* (g >= 0)));
end