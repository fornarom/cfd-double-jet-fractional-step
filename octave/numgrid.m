function G = numgrid(R,n)
% Minimal replacement of MATLAB numgrid(S,n) for GNU Octave: square grid,
% interior points numbered column-wise, boundary = 0.
G = zeros(n,n); k = 0;
for j = 2:n-1, for i = 2:n-1, k = k+1; G(i,j) = k; end, end
end
