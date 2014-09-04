%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Applied Numerical Methods
%
% Example:	13.1 Postprocessing
%
% Problem:	dphi/dt = v dphi/dx
%
% Method:	Finite Difference Method with second order central differences
%			and Fourth Order Runge-Kutta Method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

load 'Shallow_Water.data';

% Simulation parameters
 x_min           =  0.00;
 x_max           = 100.00;
 y_min           = 0.00;
 y_max           = 100.00;
 t_min           = 0.00;
 t_max           = 100.00;
     Delta_x         =  0.4;            % Experiment with 0.02, 0.05, 0.10
    Delta_y         =  0.4;            % Experiment with 0.02, 0.05, 0.10
    Delta_t         =  0.04;            % Experiment with 0.01, 0.02, 0.10

 
 %[N_t N_x]   = size(Example13_1);
  x               = x_min:Delta_x:x_max;
  y               = y_min:Delta_y:y_max;
    
 N_x             = length(x);
 N_y             = length(y);
 N_t             = (t_max-t_min)/Delta_t+1;


% Allocate arrays
%x           = linspace(x_min, x_max, N_x);
%t           = linspace(t_min, t_max, N_t);
[x y]           = meshgrid(x_min:Delta_x:x_max, y_min:Delta_y:y_max);

% Set up animated visualization of solution
figure('WindowStyle', 'docked');
colormap hot;
axes;
Solution = surf(x, y, Shallow_Water(1:N_x,:));
%axis('equal', [x_min x_max y_min y_max 1 2]);
ax=axis;
grid on;
xlabel('x');
ylabel('y');
zlabel('\phi');
view([45 25]);
    
drawnow;

for l=1:N_t-1
    rmin = l *  N_x+1;
    rmax = (l+1)*N_x;
    %Plot the solution
    set(Solution, 'ZData', Shallow_Water(rmin:rmax,:));
    axis(ax);
    title(['t = ' num2str(l+1)]);
    drawnow;        
%    pause(0.2);
end