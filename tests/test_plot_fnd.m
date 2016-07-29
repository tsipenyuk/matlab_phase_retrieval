a = -5:1:5;
b = a;
[A,B] = meshgrid(a,b);
Z = complex(A,B);

% Complex differentiable
F = Z.*Z;
G = 2*Z;
plot_fnd(Z,F,G,'z^2','2z');

% Wirtinger
ax = subplot(2,2,1);
F = abs(Z).^2;
G = conj(Z);
plot_fnd(Z,F,G,'|z|^2', '$d/dz |z|^2', ax);
title('d/dz')

ax = subplot(2,2,2);
GC = Z;
plot_fnd(Z,F,GC,'|z|^2', 'd/dz |z|^$', ax);
title('d/d conj(z)')

ax = subplot(2,2,3);
GS = 2*GC;
plot_fnd(Z,F,GS,'|z|^2', 'd/d {conj(z)} |z|^2', ax);
title('d/d real(z) + 1j * d/d imag(z)')

ax = subplot(2,2,4);
G = 2*real(Z);
plot_fnd(Z,F,G,'|z|^2', 'Fienup d/dz', ax);
title('Fienup d/dz')
