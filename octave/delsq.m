function D = delsq(G)
% Minimal replacement of MATLAB delsq for GNU Octave: 5-point discrete
% Laplacian (4 on the diagonal, -1 for each neighbour) on the grid G.
N = max(G(:)); [m,n] = size(G); I=[];J=[];S=[];
for j=1:n, for i=1:m
  k=G(i,j); if k==0, continue; end
  I(end+1)=k; J(end+1)=k; S(end+1)=4;
  nb=[i-1 j; i+1 j; i j-1; i j+1];
  for q=1:4, kk=G(nb(q,1),nb(q,2)); if kk>0, I(end+1)=k; J(end+1)=kk; S(end+1)=-1; end, end
end, end
D = sparse(I,J,S,N,N);
end
