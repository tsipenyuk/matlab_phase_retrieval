function ut = laxFriedrichs2d(u0, v0)


% $f = u \cdot v$
% f1 = u * v1
v1 = v(:,:,1);
v2 = v(:,:,2);

f1 = u .* v1;
f2 = u .* v2;
