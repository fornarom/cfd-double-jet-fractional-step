%% Double jet flow in a doubly periodic domain
% Fractional Step Procedure (projection method) for the 2D incompressible
% Navier-Stokes equations, Harlow-Welch staggered grid, 4-stage Runge-Kutta
% time integration.
%
% Author: Michele Fornaro - Computational Fluid Dynamics, MSc in Mathematical
% Engineering, University of Naples Federico II (A.Y. 2022-2023)

clear; close all; clc;

global Nx Ny Re h hsq                   % Global declaration of variables
warning off                             % Suppress warning messages

%% Domain and mesh: [0, 2*pi] x [0, 2*pi]
Lx  = 2*pi;   AR = 1;   Ly = AR*Lx;     % Domain extension
Ncx = 180;    Ncy = AR*Ncx;             % Number of discretization cells
Nx  = Ncx + 1;   Ny = Ncy + 1;          % Number of discretization points
x   = linspace(0,Lx,Nx);                % Mesh on the x direction
y   = linspace(0,Ly,Ny);                % Mesh on the y direction
hx  = x(2) - x(1);                      % Mesh size on x (uniform)
hy  = y(2) - y(1);                      % Mesh size on y (uniform)

% Uniform and homogeneous spatial mesh on a square domain
h   = hx;
hsq = h*h;    N = Nx;
Re  = Inf;                              % Inviscid flow (viscous term not used)
T   = 10;                               % Final time

% Velocity components and pressure field
U = zeros(Nx-1,Ny-1);   V = zeros(Nx-1,Ny-1);   P = zeros(Nx-1,Ny-1);

%% Initial conditions: double shear layer + small transverse perturbation
delta = pi/15;                          % Shear layer thickness
ep    = 0.05;                           % Perturbation amplitude

for j = 1:Ny-1
    if y(j) <= pi
        U(:,j) = -tanh((y(j) - pi/2)/delta);
    else
        U(:,j) = -tanh((3*pi/2 - y(j))/delta);
    end
end

for i = 1:Nx-1
    V(i,:) = -ep*sin(x(i));
end

%% Time step
Uref = 1;                               % Reference velocity
C    = 0.4;                             % Maximum Courant number
Dt   = C*h/Uref;                        % Dt for stability
Nt   = round(T/Dt);                     % Number of time steps

%% Discrete Laplacian with periodic boundary conditions
G   = numgrid('S',N+1);
Lap = -delsq(G)/hsq;

for i = 2:Nx
    j = 2;    ks = G(i,j);  kk = G(i,Ny);  Lap(ks,kk) = 1/hsq;
    j = Ny;   kn = G(i,j);  kk = G(i,2);   Lap(kn,kk) = 1/hsq;
end
for j = 2:Ny
    i = 2;    kw = G(i,j);  kk = G(Nx,j);  Lap(kw,kk) = 1/hsq;
    i = Nx;   ke = G(i,j);  kk = G(2,j);   Lap(ke,kk) = 1/hsq;
end

%% Unsteady time integration: Fractional Step + RK4
for it = 1:Nt
    t = it*Dt;

    % STAGE 1
    U1 = U;   V1 = V;
    [Fu1,Fv1] = RHS_Periodic_HW(U1,V1);

    % STAGE 2: starred velocity, pressure Poisson equation, projection
    U2  = U + Dt*0.5*Fu1;
    V2  = V + Dt*0.5*Fv1;
    DIV = Div_Calc_periodic_HW(U2,V2);
    P   = Lap\DIV(:);
    P   = reshape(P,Nx-1,Ny-1);
    [Px,Py] = Grad_Calc_periodic_HW(P);
    U2  = U2 - Px;
    V2  = V2 - Py;
    [Fu2,Fv2] = RHS_Periodic_HW(U2,V2);

    % STAGE 3
    U3  = U + Dt*0.5*Fu2;
    V3  = V + Dt*0.5*Fv2;
    DIV = Div_Calc_periodic_HW(U3,V3);
    P   = Lap\DIV(:);
    P   = reshape(P,Nx-1,Ny-1);
    [Px,Py] = Grad_Calc_periodic_HW(P);
    U3  = U3 - Px;   V3 = V3 - Py;
    [Fu3,Fv3] = RHS_Periodic_HW(U3,V3);

    % STAGE 4
    U4  = U + Dt*Fu3;
    V4  = V + Dt*Fv3;
    DIV = Div_Calc_periodic_HW(U4,V4);
    P   = Lap\DIV(:);
    P   = reshape(P,Nx-1,Ny-1);
    [Px,Py] = Grad_Calc_periodic_HW(P);
    U4  = U4 - Px;
    V4  = V4 - Py;
    [Fu4,Fv4] = RHS_Periodic_HW(U4,V4);

    % FINAL CORRECTION STEP
    U   = U + Dt*((1/6)*(Fu1 + Fu4) + (1/3)*(Fu2 + Fu3));
    V   = V + Dt*((1/6)*(Fv1 + Fv4) + (1/3)*(Fv2 + Fv3));
    DIV = Div_Calc_periodic_HW(U,V);
    P   = Lap\DIV(:);
    P   = reshape(P,Nx-1,Ny-1);
    [Px,Py] = Grad_Calc_periodic_HW(P);
    U   = U - Px;   V = V - Py;

    %% Output every 10 time steps
    if mod(it,10) == 0

        % Mass conservation check
        MaxDiv = max(max(abs(Div_Calc_periodic_HW(U,V))));
        disp(['Maximum value of Div(V) = ',num2str(MaxDiv)])

        % Total kinetic energy
        tot_en = Kinetic_energy(U,V);
        disp(['Value of the total Kinetic energy = ',num2str(tot_en)])

        tot_dv_vec(it,1) = MaxDiv;
        tot_en_vec(it,1) = tot_en;
        t_vec(it,1)      = t;

        % Interpolation at cell corners of U along y
        Iu(1:Nx-1,2:Ny-1) = (U(1:Nx-1,2:Nx-1) + U(1:Nx-1,(2:Nx-1)-1))/2;
        Iu(1:Nx-1,1)      = (U(1:Nx-1,Nx-1) + U(1:Nx-1,1))/2;
        Iu(1:Nx-1,Nx)     = Iu(1:Nx-1,1);
        Iu(Nx,1:Ny)       = Iu(1,1:Nx);

        % Interpolation at cell corners of V along x
        Iv(2:Ny-1,1:Nx-1) = (V(2:Nx-1,1:Nx-1) + V((2:Nx-1)-1,1:Nx-1))/2;
        Iv(1,1:Nx-1)      = (V(Nx-1,1:Nx-1) + V(1,1:Nx-1))/2;
        Iv(Nx,1:Nx-1)     = Iv(1,1:Nx-1);
        Iv(1:Ny,Nx)       = Iv(1:Nx,1);

        % 2D mesh
        X = repmat(x,Ny,1);   Y = repmat(y',1,Nx);

        % Velocity field
        figure(1); clf;
        if exist('streamslice','file')
            streamslice(X,Y,Iu',Iv');
        else                            % GNU Octave fallback
            quiver(X(1:6:end,1:6:end),Y(1:6:end,1:6:end), ...
                   Iu(1:6:end,1:6:end)',Iv(1:6:end,1:6:end)');
        end
        axis([0 Lx 0 Ly]); axis square;
        title(['Doubly periodic jet flow. t = ',num2str(t)]);
        drawnow;

        % Streamlines from the stream function
        figure(2); clf;
        [PSI,ZITA] = GivePsiZita_periodic(U,V,Lap);
        vneg = linspace(min(min(PSI)),0,10);
        vpos = linspace(0,max(max(PSI)),10);
        contour(x(1:Ncx),y(1:Ncx),PSI',vneg,'k'); axis square; hold on;
        contour(x(1:Ncx),y(1:Ncx),PSI',vpos,'r');
        title('Harlow-Welch. Doubly periodic jet flow. Streamlines from \psi')

        % Isovorticity contours
        figure(3); clf;
        contour(y(1:Ncx),x(1:Ncx),ZITA',[-6 -5 -4 -3 -2 -1 1 2 3 4 5 6],'k');
        axis square; hold on;
        title(['Doubly periodic jet flow. Isovorticity contours \zeta, t = ',num2str(t)]);
        drawnow; hold off

        % Enstrophy
        Z       = (ZITA.*ZITA)*h^2;
        tot_ens = sum(Z(:))/2;
        disp(['Value of the total enstrophy = ',num2str(tot_ens)])
        tot_ens_vec(it,1) = tot_ens;
    end
end

%% Conservation diagrams
% Keep only the time steps where outputs were stored
idx         = find(t_vec ~= 0);
t_vec       = t_vec(idx);
tot_en_vec  = tot_en_vec(idx);
tot_ens_vec = tot_ens_vec(idx);
tot_dv_vec  = tot_dv_vec(idx);

% Total kinetic energy error
figure(4)
plot(t_vec, tot_en_vec - tot_en_vec(1));
xlabel('Time'); ylabel('Kinetic Energy Error');

% Total enstrophy
figure(5)
plot(t_vec, tot_ens_vec);
xlabel('Time'); ylabel('Total Enstrophy');

% Maximum divergence of velocity
figure(6)
plot(t_vec, tot_dv_vec);
xlabel('Time'); ylabel('Divergence of Velocity');
