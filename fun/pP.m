function res = pP(g)
% POSITIVITY_PROJECTION Replace negative entries of g by 0
    res = real(abs(g .* (g >= 0)));
end