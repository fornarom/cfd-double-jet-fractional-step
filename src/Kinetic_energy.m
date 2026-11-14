function tot_en = Kinetic_energy(U,V)
% Total kinetic energy of the flow field (velocities interpolated at
% cell centers).

global Nx hsq

% Interpolation at cell centers
Iu(1:Nx-2,1:Nx-1) = (U(2:Nx-1,1:Nx-1) + U(1:Nx-2,1:Nx-1))/2;
Iu(Nx-1,1:Nx-1)   = (U(1,1:Nx-1) + U(Nx-1,1:Nx-1))/2;

Iv(1:Nx-1,1:Nx-2) = (V(1:Nx-1,2:Nx-1) + V(1:Nx-1,1:Nx-2))/2;
Iv(1:Nx-1,Nx-1)   = (V(1:Nx-1,1) + V(1:Nx-1,Nx-1))/2;

E   = Iu + Iv;
Esq = E.*E;
tot_en = sum(Esq(:))*hsq/2;

end
