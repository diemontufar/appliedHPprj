%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Applied Numerical Methods
%
%
% Problem:		dphi_1/dt = lambdaim*phi  
%
% Method:		Fourth Order Runge-Kutta Method Error Analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Example10_1()

    clear all; close all; clc;

    % Simulation parameters
    t_min       =  0.0;
    t_max       = 50.0;
    Delta_t     =  0.05;     % Experiment with 0.01, 0.05, 0.1
    
    phi_min     = exp0;
    N_t         = (t_max-t_min)/Delta_t+1;

    % Allocate arrays
    t           = zeros(1,   N_t);
    phi         = zeros(N_t);
    phireal = zeros(N_t);

    % Set initial condition
    t(1)        = t_min;
    phi(1)	= phi_min;

    % Time marching loop
    for l = 1:N_t-1

        t(l+1)      = t(l) + Delta_t;

        k1			= f(phi(l),                   t(l)          );
        k2			= f(phi(l)+Delta_t/2.*k1, t(l)+Delta_t/2);
        k3			= f(phi(l)+Delta_t/2.*k2, t(l)+Delta_t/2);
        k4			= f(phi(l)+Delta_t   *k3, t(l)+Delta_t  );

        phi(l+1)  = phi(l) + Delta_t*(k1/6 + k2/3 + k3/3 + k4/6);
        phireal(l) = exp(i*l);
    end

    % Plot the solution
    figure('WindowStyle', 'docked');
    plot(t, phi,  t , [exp(i*t)]);
    %plot(t, [exp(i*t)], '.-');
    axis([t_min t_max -1.5 1.5]);
    grid on;
    xlabel('t');
    ylabel('\phi');
    legend('\phi_1','\phireal','location','SouthEast');

    figure('WindowStyle', 'docked');
%    plot(t, phi,  t , [exp(i*t)]);
    plot(t, [exp(i*t)], '.-');
    axis([t_min t_max -1.5 1.5]);
    grid on;
    xlabel('t');
    ylabel('\phi');
    legend('\phireal','location','SouthEast');

    figure('WindowStyle', 'docked');
%    plot(t, phi,  t , [exp(i*t)]);
    plot(t, phi, '.-');
    axis([t_min t_max -1.5 1.5]);
    grid on;
    xlabel('t');
    ylabel('\phi');
    legend('\phi','location','SouthEast');

return

function k = f(phi, t)

    k     = i*phi;

return

    





