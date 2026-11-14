function [PSI, ZITA] = GivePsiZita_periodic(U,V,Lap)
% Vorticity (zeta) and stream function (psi) from the velocity field,
% doubly periodic domain.

global h Nx Ny

% Interpolation at cell centers
Iu(1:Nx-2,1:Nx-1) = (U(2:Nx-1,1:Nx-1) + U(1:Nx-2,1:Nx-1))/2;
Iu(Nx-1,1:Nx-1)   = (U(1,1:Nx-1) + U(Nx-1,1:Nx-1))/2;

Iv(1:Nx-1,1:Nx-2) = (V(1:Nx-1,2:Nx-1) + V(1:Nx-1,1:Nx-2))/2;
Iv(1:Nx-1,Nx-1)   = (V(1:Nx-1,1) + V(1:Nx-1,Nx-1))/2;

% Velocity derivatives
Iuy(1:Nx-1,1:Ny-2) = (Iu(1:Nx-1,2:Nx-1) - Iu(1:Nx-1,(2:Nx-1)-1))/h;
Iuy(1:Nx-1,Ny-1)   = (Iu(1:Nx-1,1) - Iu(1:Nx-1,Nx-1))/h;

Ivx(1:Ny-2,1:Nx-1) = (Iv(2:Nx-1,1:Nx-1) - Iv((2:Nx-1)-1,1:Nx-1))/h;
Ivx(Ny-1,1:Nx-1)   = (Iv(1,1:Nx-1) - Iv(Nx-1,1:Nx-1))/h;

% Vorticity field
ZITA = -Iuy + Ivx;

% Stream function from the Poisson equation
Psi = Lap\ZITA(:);
PSI = reshape(Psi,Nx-1,Ny-1);

end
