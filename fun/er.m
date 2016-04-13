function [g_new, error] = er(g, A)
% ER Error reduction algorithm
% <TEST THIS>
    g_new = pP(pM(g, A));
    error = eM(g_new, A);
end