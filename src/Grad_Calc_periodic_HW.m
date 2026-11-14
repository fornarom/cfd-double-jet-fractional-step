function [Px,Py] = Grad_Calc_periodic_HW(P)
% Gradient of the pressure field at velocity nodes, Harlow-Welch staggered
% grid, doubly periodic domain.

global h Nx Ny

Px = zeros(Nx-1,Ny-1);   Py = zeros(Nx-1,Ny-1);

% Internal points
Px(2:Nx-1,1:Nx-1) = (P(2:Nx-1,1:Nx-1) - P((2:Nx-1)-1,1:Nx-1))/h;
Py(1:Nx-1,2:Nx-1) = (P(1:Nx-1,2:Nx-1) - P(1:Nx-1,(2:Nx-1)-1))/h;

% West boundary (Px) and North boundary (Py)
Px(1,:)      = (P(1,:) - P(Nx-1,:))/h;
Py(1:Nx-1,1) = (P(:,1) - P(:,Nx-1))/h;

end
