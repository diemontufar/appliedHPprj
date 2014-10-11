%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Applied Numerical Methods
%
% Assignment 1:	The Shallow water Equation
% 
% Authors: Diego Montufar / Andres Chaves
%
% Date: 26/09/2014
%
% Method:		Stability and error Analysis of 4 order Runge-Kutta methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Stability_4RK()

   clear all; close all; clc;
    
    global v Delta_x N_x;

    % Simulation parameters
    x_min           =  0.00;
    x_max           = 10.00;
    Delta_x         =  0.5;            % Experiment with 0.02, 0.05, 0.10
    Delta_t         =  0.01;           % Experiment with 0.01, 0.02, 0.10
    x               =  x_min:Delta_x:x_max;
    N_x             =  length(x);
    
    % Allocate arrays
    [X, Y]          = meshgrid(-4:0.1:4, -4:0.1:4);
    Z               = X + i*Y;
    sigma           = abs(1 + Z + (Z.^2)/2 + (Z.^3)/6 + (Z.^4)/24);
    
    % Plot the stability diagram
     figure('WindowStyle', 'docked');
     %[Xi Lambda]     = eig(-v/(2*Delta_x).*(-1.*diag(ones(N_x-1, 1), -1) + diag(ones(N_x-1, 1), 1)));
     contourf(X, Y, sigma, [1 1],'-k');
     hold on;
     %plot(real(diag(sigma)*Delta_t),imag(diag(sigma))*Delta_t, '.', 'MarkerSize', 20);
     axis('equal', [-4 4 -4 4]);
     grid on;
     xlabel('\lambda_{Re}\Delta t');
     ylabel('\lambda_{Im}\Delta t');

return
    


