function DIV = Div_Calc_periodic_HW(U,V)
% Divergence of the velocity field at pressure nodes (cell centers),
% Harlow-Welch staggered grid, doubly periodic domain.

global h Nx

DIV = nan(Nx-1,Nx-1);

% Internal points
DIV(1:Nx-2,1:Nx-2) = (U((1:Nx-2)+1,1:Nx-2) - U(1:Nx-2,1:Nx-2) + ...
                      V(1:Nx-2,(1:Nx-2)+1) - V(1:Nx-2,1:Nx-2))/h;
% South boundary
DIV(1:Nx-2,Nx-1) = (U((1:Nx-2)+1,Nx-1) - U(1:Nx-2,Nx-1) + ...
                    V(1:Nx-2,1) - V(1:Nx-2,Nx-1))/h;
% East boundary
DIV(Nx-1,1:Nx-2) = (U(1,1:Nx-2) - U(Nx-1,1:Nx-2) + ...
                    V(Nx-1,(1:Nx-2)+1) - V(Nx-1,1:Nx-2))/h;
% South-East corner
DIV(Nx-1,Nx-1) = (U(1,Nx-1) - U(Nx-1,Nx-1) + V(Nx-1,1) - V(Nx-1,Nx-1))/h;

end
