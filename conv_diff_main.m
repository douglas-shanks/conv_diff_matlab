clc
clear all
% main driver for convection diffusion problem

% grid points
N = 50;
% grid spacing
dx = 1./(N-1);
dx2 = dx*dx;
% final time
T = 0.04;


% initial solution
[X,Y] = meshgrid([0:dx:1],[0:dx:1]);
u0 = sin(X).*sin(Y);
figure(1)
surf([0:dx:1],[0:dx:1],u0)
title('Initial solution')
xlabel('x')
ylabel('y')
colorbar
shading interp
axis tight
view(2)
% reshape initial solution for algorithm
u1 = reshape(u0,N*N,1);

%% Explicit method
% number of time steps
Nt = 1000;
% time step
dt = T/Nt;

% convection coefficient
alpha = 1.0;
Rc = alpha/(2.0*(dx + dx));
% diffusion coefficient
beta = 1.0;
Rd = beta*dt./((dx2 + dx2));
% Peclet number
P = alpha*(dx+dx)/(2.0*beta);
% Stability check for explici
if (Rd>=0.25 || P >=1 )
  'Explicit method unstable'
endif

% explicit matrix
A_exp = conv_diff_exp(N,Rd,P);

% main loop
for i = 2:Nt
  u = A_exp*u1;
  u1 = u;
end
figure(2)
surf([0:dx:1],[0:dx:1],reshape(u,N,N))
title(sprintf('Solution, with explicit method, at time %0.5e',T))
xlabel('x')
ylabel('y')
colorbar
shading interp
axis tight
view(2)

% Implicit for diffusion, explicit convection method
Nt = 10
dt = T/Nt;
% convection coefficient
alpha = 1.0;
Rc = alpha*dt/(2.0*(dx + dx));
% diffusion coefficient
beta = 1.0;
Rd = beta*dt./((dx2 + dx2));
% Peclet number
P = alpha*(dx+dx)/(2.0*beta);

% implicit matrix
A_imp = conv_diff_imp(N,Rd);
A_rhs_imp = rhs_imp(N,Rc);
b_imp = zeros(N*N,1);
% initial solution
[X,Y] = meshgrid([0:dx:1],[0:dx:1]);
u0 = sin(X).*sin(Y);

% reshape initial solution for algorithm
u1 = reshape(u0,N*N,1);

% main loop
for i = 2:Nt
  b_imp = A_rhs_imp*u1;
  u = pcg(A_imp,b_imp,1e-6,1000);
  u1 = u;
end

figure(3)
surf([0:dx:1],[0:dx:1],reshape(u,N,N))
title(sprintf('Solution, with implicit method, at time %0.5e',T))
xlabel('x')
ylabel('y')
colorbar
shading interp
axis tight
view(2)