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
% Method:		Finite Difference Method with fourth order central differences
%			    and Fourth Order Runge-Kutta Method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Shallow_Water()

    clear all; close all; clc;
    
    global g Delta_x Delta_y N_x N_y;

    % Simulation parameters
    x_min           =  0.00;
    x_max           = 100.00;
    y_min           = 0.00;
    y_max           = 100.00;
    t_min           = 0.00;
    t_max           = 100.00;
    Delta_x         =  0.5;           
    Delta_y         =  0.5;            
    Delta_t         =  0.1;            
    
     x               = x_min:Delta_x:x_max;
     y               = y_min:Delta_y:y_max;
     v               =  1.00;
    
    N_x             = length(x);
    N_y             = length(y);
    N_t             = (t_max-t_min)/Delta_t+1;
    
    % Allocate arrays
    phiVx            = zeros(N_x, N_y);
    phiVy            = zeros(N_x, N_y);
    phiH             = zeros(N_x, N_y);
    [x y]            = meshgrid(x_min:Delta_x:x_max, y_min:Delta_y:y_max);
        
    % Set initial condition
    phiH(:,:) = 1 + 0.5*exp(-1/25*((x-30).^2+(y-30).^2));
    g = 9.81;
    
    % Set up animated visualization of solution
    figure('WindowStyle', 'docked');
    colormap winter;
    axes;
    Solution = surf(x, y, phiH(:,:));
    ax=axis;
    grid on;
    xlabel('x');
    ylabel('y');
    zlabel('phiH');
    view([45 25]);
    
    drawnow;

    % Time marching loop
    for l=1:N_t-1

        [k1Vx, k1Vy, k1h]  = f(phiVx(:,:), phiVy(:,:), phiH(:,:));
        [k2Vx, k2Vy, k2h]  = f(phiVx(:,:) + Delta_t/2*k1Vx, phiVy(:,:) + Delta_t/2*k1Vy, phiH(:,:) + Delta_t/2*k1h);
        [k3Vx, k3Vy, k3h]  = f(phiVx(:,:) + Delta_t/2*k2Vx, phiVy(:,:) + Delta_t/2*k2Vy, phiH(:,:) + Delta_t/2*k2h);
        [k4Vx, k4Vy, k4h]  = f(phiVx(:,:) + Delta_t*k3Vx, phiVy(:,:) + Delta_t*k3Vy, phiH(:,:) + Delta_t*k3h);
        
        phiVx(:,:)  =   phiVx(:,:) + Delta_t  *(k1Vx/6 + k2Vx/3 + k3Vx/3 + k4Vx/6);
        phiVy(:,:)  =   phiVy(:,:) + Delta_t  *(k1Vy/6 + k2Vy/3 + k3Vy/3 + k4Vy/6);
        phiH(:,:)   =   phiH(:,:) + Delta_t  *(k1h/6 + k2h/3 + k3h/3 + k4h/6);

       % Plot the solution
       set(Solution, 'ZData', phiH(:,:));
       axis(ax);
       title(['t = ' num2str(l+1)]);
       drawnow;
       
   
    end

    return

function [kVx,kVy,kH] = f(phiVx,phiVy,phiH)

    global g Delta_x Delta_y N_x N_y;
    
    kVx = zeros(N_x,N_y);
    kVy = zeros(N_x,N_y);
    kH = zeros(N_x,N_y);
    
    for i=1:N_x;
        for j=1:N_y;
            
            % Periodic Boundary condition
             ip1=i+1;
             ip2=i-1;
             ip3=i+2;
             ip4=i-2;
             jp1=j+1;
             jp2=j-1;
             jp3=j+2;
             jp4=j-2;
            
            if(i==N_x)
                 ip1=1;
                 ip3=2;
            end
            
            if(i==N_x-1)
                 ip3=1;
            end
            
            if(i==1)
                ip2=N_x;
                ip4=N_x-1;
            end
            
            if(i==2)
                ip4=N_x;
            end
            
            if(j==N_y)
                  jp1=1;
                  jp3=2;
            end
            
            if(j==N_y-1)
                 jp3=1;
            end
            
            if(j==1)
                jp2=N_y;
                jp4=N_y-1;
            end
            
            if(j==2)
                jp4=N_y;
            end
            
            % Function evaluation 4th order central differences
            kVx(i,j)= (-g/(12*Delta_x)*(phiH(ip4,j)-8*(phiH(ip2,j))+8*phiH(ip1,j)-phiH(ip3,j))) - (phiVx(i,j)/(12*Delta_x)*(phiVx(ip4,j)-8*(phiVx(ip2,j))+8*phiVx(ip1,j)-phiVx(ip3,j))) - (phiVx(i,j)/(12*Delta_y)*(phiVx(i,jp4)-8*(phiVx(i,jp2))+8*phiVx(i,jp1)-phiVx(i,jp3)));
            kVy(i,j)= (-g/(12*Delta_y)*(phiH(i,jp4)-8*(phiH(i,jp2))+8*phiH(i,jp1)-phiH(i,jp3))) - (phiVx(i,j)/(12*Delta_x)*(phiVy(ip4,j)-8*(phiVy(ip2,j))+8*phiVy(ip1,j)-phiVy(ip3,j))) - (phiVy(i,j)/(12*Delta_y)*(phiVy(i,jp4)-8*(phiVy(i,jp2))+8*phiVy(i,jp1)-phiVy(i,jp3)));
            kH(i,j) = (-phiVx(i,j)/(12*Delta_x)*(phiH(ip4,j)-8*(phiH(ip2,j))+8*phiH(ip1,j)-phiH(ip3,j))) - (phiH(i,j)/(12*Delta_x)*(phiVx(ip4,j)-8*(phiVx(ip2,j))+8*phiVx(ip1,j)-phiVx(ip3,j))) - (phiVy(i,j)/(12*Delta_y)*(phiH(i,jp4)-8*(phiH(i,jp2))+8*phiH(i,jp1)-phiH(i,jp3))) - (phiH(i,j)/(12*Delta_y)*(phiVy(i,jp4)-8*(phiVy(i,jp2))+8*phiVy(i,jp1)-phiVy(i,jp3)));
  
        end
    end

return
    


