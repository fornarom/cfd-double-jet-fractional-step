function [Fu,Fv] = RHS_Periodic_HW(U,V)
% Right Hand Side of a 2D Harlow-Welch discretization of the incompressible
% Navier-Stokes equations (convective terms, inviscid case) with periodic
% boundary conditions, for use in a projection method.

global Nx Ny h

Fu = zeros(size(U));   Fv = zeros(size(V));

%--------------------------------------------------------------------------
% U component
%--------------------------------------------------------------------------

% Term d(I(u)_x * I(u)_x)/dx
% Interpolation at cell centers of U along x
Iu(1:Nx-2,1:Ny-1) = (U((1:Nx-2)+1,1:Ny-1) + U(1:Nx-2,1:Ny-1))/2;
Iu(Nx-1,:)        = (U(1,:) + U(Nx-1,:))/2;
Iusq = Iu.*Iu;

% Difference along x at vertical cell faces
i = 2:Nx-1;   j = 1:Ny-1;
Fu(i,j)      = Fu(i,j) + (Iusq(i,j) - Iusq(i-1,j))/h;
Fu(1,1:Ny-1) = Fu(1,1:Ny-1) + (Iusq(1,1:Ny-1) - Iusq(Nx-1,1:Ny-1))/h;

% Term d(I(u)_y * I(v)_x)/dy
Iu = zeros(Nx-1,Ny-1);   Iv = zeros(Nx-1,Ny-1);

% Interpolation at cell corners of u along y
Iu(1:Nx-1,2:Nx-1) = (U(1:Nx-1,1:Nx-2) + U(1:Nx-1,(1:Nx-2)+1))/2;
Iu(1:Nx-1,1)      = (U(1:Nx-1,1) + U(1:Nx-1,Nx-1))/2;

% Interpolation at cell corners of v along x
Iv(2:Nx-1,1:Nx-1) = (V(2:Nx-1,1:Nx-1) + V((2:Nx-1)-1,1:Nx-1))/2;
Iv(1,1:Nx-1)      = (V(1,1:Nx-1) + V(Nx-1,1:Nx-1))/2;

IuIv = Iu.*Iv;

Fu(1:Nx-1,1:Ny-2) = Fu(1:Nx-1,1:Ny-2) + ...
    (IuIv(1:Nx-1,(1:Ny-2)+1) - IuIv(1:Nx-1,1:Ny-2))/h;
Fu(1:Nx-1,Nx-1)   = Fu(1:Nx-1,Nx-1) + ...
    (IuIv(1:Nx-1,1) - IuIv(1:Nx-1,Nx-1))/h;

%--------------------------------------------------------------------------
% V component
%--------------------------------------------------------------------------

% Term d(I(v)_y * I(v)_y)/dy
Iv = zeros(size(V));
Iv(:,1:Ny-2) = (V(:,(1:Ny-2)+1) + V(:,1:Ny-2))/2;
Iv(:,Ny-1)   = (V(:,Ny-1) + V(:,1))/2;
Ivsq = Iv.*Iv;

Fv(1:Nx-1,2:Ny-1) = Fv(1:Nx-1,2:Ny-1) + ...
    (Ivsq(1:Nx-1,2:Ny-1) - Ivsq(1:Nx-1,(2:Ny-1)-1))/h;
Fv(1:Nx-1,1)      = Fv(1:Nx-1,1) + (Ivsq(1:Nx-1,1) - Ivsq(1:Nx-1,Ny-1))/h;

% Term d(I(u)_y * I(v)_x)/dx
Fv(1:Nx-2,1:Ny-1) = Fv(1:Nx-2,1:Ny-1) + ...
    (IuIv((1:Nx-2)+1,1:Ny-1) - IuIv(1:Nx-2,1:Ny-1))/h;
Fv(Nx-1,1:Ny-1)   = Fv(Nx-1,1:Ny-1) + (IuIv(1,1:Ny-1) - IuIv(Nx-1,1:Ny-1))/h;

end
