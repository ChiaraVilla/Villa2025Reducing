%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  "Reducing phenotype-structured PDE models of cancer evolution to   %%%
%%%       systems of ODEs: a generalised moment dynamics approach"      %%%
%%%                                                                     %%%
%%%   C. Villa, P. Maini, A. Browning, A. Jenner, S. Hamis, T. Cassidy  %%%
%%%                         Preprint (2025)                             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%% 
%%% Villa2025Reducing:                                                  %%%
%%% Code to solve the PDE model, and to construct and solve the ODE     %%%
%%% system approximating the dynamics of the moments of the phenotypic  %%%
%%% density distribution for N=2,3 and M=1,2,3 and to compare the       %%%
%%% solutions                                                           %%%
%%%                                                                     %%%
%%% Copyright:  C. Villa (*)                                            %%%
%%%                                                                     %%%
%%% (*) Email: chiara.villa[at]inria.fr                                 %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%% 
%%% Note 1:                                                             %%% 
%%% the code is set up to test any f and V: simply change the           %%% 
%%% functions V_DEF(x,par) and f_DEF(x,par). If your functions are      %%%
%%% time dependent you many need to edit these to allow for the extra   %%%
%%% input variable, along with functions dnf_dxn(M,par) and             %%%
%%% dnV_dxn(M,par), and fn =@(n,mi) and Vn =@(n,mi) in the code.        %%%
%%%                                                                     %%% 
%%% Note 2:                                                             %%% 
%%% The code uses explicit schemes to solve the PDE and ODE system.     %%%
%%% This my not be the best choice for custom f and V and for all       %%%
%%% combinations of N and M, and the scheme may be unstable in certain  %%%
%%% cases.                                                              %%% 
%%%                                                                     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example 2 set up
clear all
clc
date = '100325'; 

%% Problem set up

%%% Model parameters - change directly in 'parameters.m'
par = parameters;

%%% Phenotypic domain discretisation
par.res = 1.2;   % Rescale factor
par.a = 0;          % Lower bound
par.b = 1.2/par.res;           % Upper bound
par.Nx = 1000;        % Number of grid cell volumes
par.dx = (par.b-par.a)/par.Nx; % Grid cell size

xint = linspace(par.a,par.b,par.Nx+1);     % Grid cell interfaces
x(1:par.Nx) = xint(1:par.Nx) + 0.5*par.dx; % Grid cell centers

%%% Times at which we want the solution
Tend = 35;
resol = 20; % image resolution
tspan = linspace(0,Tend,Tend*resol+1);
dt=0.001;

%%% CHOOSE APPROXIMATION PARAMETERS
N=2; % Number of moments chosen to track
M=2; % Number of derivatives in Taylor series approximation
if M>3 || N>3
    error('This code is set up only for N<4 and M<4')
elseif N==1
    error('This code is not set up for N=1 (cannnot test gaussian closure)')
end

%%% CHOOSE SYSTEM CLOSURE
closure = 'gaussian';
%closure = 'truncate';
if closure~='gaussian' & closure~='truncation'
    error('Choose between gaussian or truncation closure')
end

%%% Initial condition
p = exp(-((x-(par.x0/par.res))./(par.sig/par.res)).^2);
mk = Moments(p,x,par.dx,N);
P = mk(1);

%% Solve PDE and store data

%%% Initialise functions and storing variables
Vint = V_DEF(xint,par); 
f = f_DEF(x,par); 
DMx = (par.beta/(par.res^2))*DM_def(par.Nx,par.dx);
pstore = [p];
mk_PDE_store = [mk];

%%% Iterate in time (explicit scheme)
for i=1:1:Tend/dt
    adv = Upwind(p,Vint,par.Nx,par.dx);
    diff=(DMx*p')';
    reac = (f - P/par.k).*p;
    p = p + dt*( diff - adv + reac );
    P = sum(p)*par.dx;
    if mod(i,(tspan(2)-tspan(1))/dt)==0 %any(tspan==i*dt) 
        i
        plot(x,p)
        drawnow
        pstore=[pstore;p];
        mk_PDE_store = [mk_PDE_store;Moments(p,x,par.dx,N)];
    end 
end
save(strcat('Saved_PDE_sol',date,'.mat'),'pstore','mk_PDE_store','x','tspan');

%% Solve ODE system
 
%%%APPROXIMATE f and V (here up to M=3)
dVn_fun = dnV_dxn(3,par); 
dfn_fun = dnf_dxn(3,par);
fn =@(n,mi) (n>M)*0 + (n<=M)*dfn_fun{n+1}(mi);
Vn = @(n, mi) (n>M)*0 + (n<=M)*dVn_fun{n+1}(mi);

switch N
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 % N=2
        %%% Initial condition 
        P = mk(1);
        m1 = mk(2);
        m2 = mk(3);
        %%% Storage
        mk_ODE_store = [mk];
        %%% Initalise closure
        if closure=='gaussian' 
            m3 = m1*m2 + 3*(m2-m1^2)*m1;
            m4 = m1*m3 + 4*(m2-m1^2)*m2;
            m5 = m1*m4 + 5*(m2-m1^2)*m3;
        else % truncation closure
            m3 = 0;
            m4 = 0;
            m5 = 0;
        end
        %%% Iterate in time (explicit scheme)
        for i=1:1:Tend/dt
            %%% ODE system recurrent coefficients
            C1 = m2-m1^2; 
            C2 = m3-3*m1*m2+2*m1^3;
            C3 = m3-2*m1*m2+m1^3;
            C4 = m4-3*m1*m3+3*(m1^2)*m2-m1^4;
            C5 = m3-m1*m2;
            C6 = m4-2*m1*m3+(m1^2)*m2;
            C7 = m5-3*m1*m4+3*(m1^2)*m3-(m1^3)*m2;
            %%% ODE system RHS
            dPdt = ( fn(0,m1) +fn(2,m1)*C1 +fn(3,m1)*C2 )*P -(P^2)/par.k;
            dm1dt = Vn(0,m1) + (fn(1,m1)+Vn(2,m1)-fn(2,m1)*m1)*C1 +fn(3,m1)*C4;
            dm2dt = -m2*( fn(0,m1) +fn(2,m1)*C1 +fn(3,m1)*C2 ) +2*par.beta/(par.res^2) ...
                    +2*(Vn(0,m1)*m1 +Vn(1,m1)*C1 +Vn(2,m1)*C3 +Vn(3,m1)*C4 ) ...
                    +fn(0,m1)*m2 +fn(1,m1)*C5 +fn(2,m1)*C6 +fn(3,m1)*C7;
            %%% Update moments 
            P = P + dt*dPdt;
            m1 = m1 + dt*dm1dt;
            m2 = m2 + dt*dm2dt;
            if closure=='gaussian' % If gaussian closure update HOMs
                m3 = m1*m2 + 3*(m2-m1^2)*m1;
                m4 = m1*m3 + 4*(m2-m1^2)*m2;
                m5 = m1*m4 + 5*(m2-m1^2)*m3;
            end
            %%% Store
            if mod(i,(tspan(2)-tspan(1))/dt)==0 
                mk_ODE_store = [mk_ODE_store;[P,m1,m2]];
            end 
        end
        %%% Save
        save(strcat('Saved_ODE_sol_N',num2str(N),'_M',num2str(M),'_',closure,'_',date,'.mat'),'mk_ODE_store','N','M','closure','tspan');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3 % N=3
        %%% Initial condition 
        P = mk(1);
        m1 = mk(2);
        m2 = mk(3);
        m3 = mk(4);
        %%% Storage
        mk_ODE_store = [mk];
        %%% Initalise closure
        if closure=='gaussian' 
            m4 = m1*m3 + 4*(m2-m1^2)*m2;
            m5 = m1*m4 + 5*(m2-m1^2)*m3;
            m6 = m1*m5 + 6*(m2-m1^2)*m4;
        else % truncation closure
            m4 = 0;
            m5 = 0;
            m6 = 0;
        end
        %%% Iterate in time (explicit scheme)
        for i=1:1:Tend/dt
            %%% ODE system recurrent coefficients
            C1 = m2-m1^2; 
            C2 = m3-3*m1*m2+2*m1^3;
            C3 = m3-2*m1*m2+m1^3;
            C4 = m4-3*m1*m3+3*(m1^2)*m2-m1^4;
            C5 = m3-m1*m2;
            C6 = m4-2*m1*m3+(m1^2)*m2;
            C7 = m5-3*m1*m4+3*(m1^2)*m3-(m1^3)*m2;
            C8 = m4-m1*m3;
            C9 = m5-2*m1*m4+(m1^2)*m3;
            C10 = m6-3*m1*m5+3*(m1^2)*m4-(m1^3)*m3;
            %%% ODE system RHS
            dPdt = ( fn(0,m1) +fn(2,m1)*C1 +fn(3,m1)*C2 )*P -(P^2)/par.k;
            dm1dt = Vn(0,m1) + (fn(1,m1)+Vn(2,m1)-fn(2,m1)*m1)*C1 +fn(3,m1)*C4;
            dm2dt = -m2*( fn(0,m1) +fn(2,m1)*C1 +fn(3,m1)*C2 ) +2*par.beta/(par.res^2) ...
                    +2*(Vn(0,m1)*m1 +Vn(1,m1)*C1 +Vn(2,m1)*C3 +Vn(3,m1)*C4 ) ...
                    +fn(0,m1)*m2 +fn(1,m1)*C5 +fn(2,m1)*C6 +fn(3,m1)*C7;
            dm3dt = -m3*( fn(0,m1) +fn(2,m1)*C1 +fn(3,m1)*C2 ) ...
                    +6*(par.beta/(par.res^2))*m1 ...
                    +3*( Vn(0,m1)*m2 +Vn(1,m1)*C5 +Vn(2,m1)*C6 +Vn(3,m1)*C7 ) ...
                    +fn(0,m1)*m3 +fn(1,m1)*C8 +fn(2,m1)*C9 +fn(3,m1)*C10;
            %%% Update moments 
            P = P + dt*dPdt;
            m1 = m1 + dt*dm1dt;
            m2 = m2 + dt*dm2dt;
            m3 = m3 + dt*dm3dt;
            if closure=='gaussian' % If gaussian closure update HOMs
                m4 = m1*m3 + 4*(m2-m1^2)*m2;
                m5 = m1*m4 + 5*(m2-m1^2)*m3;
                m6 = m1*m5 + 6*(m2-m1^2)*m4;
            end
            %%% Store
            if mod(i,(tspan(2)-tspan(1))/dt)==0 
                mk_ODE_store = [mk_ODE_store;[P,m1,m2,m3]];
                subplot(1,4,1)
                plot(mk_ODE_store(:,1))
                axis square
                subplot(1,4,2)
                plot(mk_ODE_store(:,2))
                axis square
                subplot(1,4,3)
                plot(mk_ODE_store(:,3))
                axis square
                subplot(1,4,4)
                plot(mk_ODE_store(:,4))
                axis square
                drawnow
            end 
        end
        %%% Save
        save(strcat('Saved_ODE_sol_N',num2str(N),'_M',num2str(M),'_',closure,'_',date,'.mat'),'mk_ODE_store','N','M','closure','tspan');
end %% Plot

%% ANNEXED FUNCTIONS

%%% Parameters
function [par] = parameters
par.fmax = 2;  % Max prolif rate
par.xm  = 0.4; % Michaelis constant
par.Vmax = 0.5; % Max velocity
par.k = 0.1; % Carrying capacity
par.beta  = 1e-04; % Diffusion coefficient
par.om1 = 2; % Velocity skewedness parameter
par.om2 = 1; % Velocity skewedness parameter
par.x0 = 0.2; % Initial mean
par.sig = 0.02; % Initial variance
end
%%% Drift velocity
function [V] = V_DEF(x,par) 
% V = (par.Vmax/par.res)*tanh((par.res*x).^par.om).*tanh(1-(par.res*x));
V = (par.Vmax/par.res)*tanh((par.res*x).^par.om1).*tanh((1-(par.res*x)).^par.om2);
end
%%% Intrinsic growth rate
function [f] = f_DEF(x,par)
f = par.fmax*x./(par.xm/par.res+x);
end

%% Approximate V and f

function dfn_fun = dnf_dxn(M,par)
% Define symbolic V
syms x
f_sym = f_DEF(x,par);
% Compute derivatives
dfn_sym = arrayfun(@(n) diff(f_sym,x,n), 0:1:M, 'UniformOutput', false);
% Convert sybmolic function to Matlab function for evaluation
dfn_fun = cellfun(@(dfn) matlabFunction(dfn, 'Vars', [x]), dfn_sym, 'UniformOutput', false);
end

function dVn_fun = dnV_dxn(M,par)
% Define symbolic V
syms x
V_sym = V_DEF(x,par);
% Compute derivatives
dVn_sym = arrayfun(@(n) diff(V_sym,x,n), 0:1:M, 'UniformOutput', false);
% Convert sybmolic function to Matlab function for evaluation
dVn_fun = cellfun(@(dVn) matlabFunction(dVn, 'Vars', [x]), dVn_sym, 'UniformOutput', false);
end

%% Calculate moments

function mk = Moments(p,x,dx,N)
P = sum(p)*dx;    % Compute total mass
p = p./P;         % Normalise p
mk = [P];         % Add zero-th moment for completeness
for k=1:N         % Compute first N moments
    mk = [mk, sum((x.^k).*p)*dx];
end
end

%% Functions to solve PDE

%%% Advection term (first order upwind)
function [dvp_dx] = Upwind(u,vim12,nx,dx)
% INPUT:
% - u: advected quantity at cell centers 
% - vim12: drift velocity at cell interfaces 
% - nx: number of grid cells 
% - dx: grid cell width (assuming uniform grid)
% Extend u to readily apply no-flux boundary conditions
u = [0,u,0];
% Flux at grid cell interfaces: first order upwind
fim12 = u(1:nx+1).*max(0,vim12) + u(2:nx+2).*min(0,vim12); 
% d(vp)/dx at cell centers (1st order centered difference) 
dvp_dx = (fim12(2:nx+1)-fim12(1:nx))./dx; 
end

%%% Diffusion matrix (finite volume scheme with no-flux BCs)
function DM_mat = DM_def(Nx,dx)
persistent DM;  % Declare DM as a persistent variable
% Check if DM exists and has the correct size
if isempty(DM) || any(size(DM) ~= [Nx, Nx])
    %%% Coefficient matrix
    DM = -2*eye(Nx);
    DM(2:Nx,1:Nx-1) = DM(2:Nx,1:Nx-1) + eye(Nx-1);
    DM(1:Nx-1,2:Nx) = DM(1:Nx-1,2:Nx) + eye(Nx-1);
    %%% No-flux boundary conditions (from finite volume scheme)
    DM(1,1) = -1;
    DM(1,2) = 1;
    DM(Nx,Nx-1) = 1;
    DM(Nx,Nx) = -1;  
    %%% Output
    DM = DM./(dx^2);
end
DM_mat = DM;
end
