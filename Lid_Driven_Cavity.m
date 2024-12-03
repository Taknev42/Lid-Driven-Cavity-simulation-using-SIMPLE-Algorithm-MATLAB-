close all;
clc ;
clear all;

% We plan to implement the SIMPLE algo on a Lid Driven Cavity Flow.

% Resolved Problems : 1. We are not exactly following an explicit scheme
% 2. The pressure corrections of the nb of boundary points should be negle.

%% Initialise 2DGrid values

global n_x
global n_y
global dx
global dy
global w
global rho
global mu 
global P_unstag 
global u_xstag 
global v_ystag

n_x = 50 ;           % number of grid points in X-dir
n_y = 50 ;           % number of grid points in Y-dir
L_x = 1 ;           
L_y = 1 ;
dx = L_x/n_x ;        % Center grid points
dy = L_y/n_y ;
w = 1 ;               % unit depth
rho = 1000 ;          % density
kin_mu = 0.01 ;       % kinematic viscosity
mu = rho*kin_mu ;     % viscosity
alpha_v = 0.5   ;     % Under- Relaxation constant for velocities
alpha_P = 0.0001 ;       % Under- Relaxation constant for P
eps = 1e-5 ;          % Convergence criteria
max_iter = 1000; 


%% Step1: Initialise Pressure and Velocity at the centers of UnStag grid

% A useful matrix for physical position to matrix grid locations
% Here we take our matrices to have the 
% topmost point @ row = 0 and bottomost point @ row = n
% Num of rows = height ; Num of columns = width
position_helper = [["" "Go Up (0,1)@(i-1,j)" ""]; ["Go Left (-1,0)@(i,j-1)" "Physical Origin (0,0)@(i,j)" "Go Right (+1,0)@(i,j+1)"] ;["" "Go Down (0,-1)@(i+1,j)" ""]] ;

% Initialising new Pressure and velcoity at all the grids
P_unstag_new = ones(n_y,n_x) ;
u_xstag_new = 0.005*ones(n_y,n_x+1 ) ;
v_ystag_new = 0.005*ones(n_y+1, n_x) ;

% Set necessary boundary values

% At the first and last columns of xStag grid
for i = 1:n_y  % row-wise
    u_xstag_new(i,1) = 0 ;
    u_xstag_new(i,n_x+1) = 0 ;
end

% At the first and last rows of the yStag grid
for j = 1:n_x           % column-wise
    v_ystag_new(1,j) = 0 ;
    v_ystag_new(n_y+1,j) = 0 ;
end


% Initialise the Pressure and Velocity corrections
% we dont need old vs new here
P_unstag_c = zeros(n_y, n_x) ;
u_xstag_c = zeros(n_y,n_x+1 ) ;
v_ystag_c = zeros(n_y+1, n_x) ;


% Initialising old Pressure and velcoity at all the grids
P_unstag = P_unstag_new ;
u_xstag = u_xstag_new ;
v_ystag = v_ystag_new ;

% Initialising pred Pressure and velcoity at all the grids
% we dont need old vs new here
P_unstag_pred = P_unstag_new ;
u_xstag_pred = u_xstag_new ;
v_ystag_pred = v_ystag_new ;

%% Iterating till convergence

for iter = 1:max_iter      %iter-wise
    
    % Step2 : Predict velocity values using old velocities in the momentum eqns

    % Starting with the Top-left-1 corner of xStag grid
    % we plug-in values at all points except the first&last columns of u_xStag matrix
    for i = 1:n_y               % row-wise             
        for j = 2:n_x               % inner column-wise

            % computing ae coefficient at (i,j) coordinate in xStag grid 
            ae = get_ae(i,j, u_xstag, v_ystag) ;     % according to old values

            P_E = P_unstag(i,j) ;                    % according to old values
            P_P = P_unstag(i,j-1) ;
            A_E = dy*w ;
            A_P = dy*w ;
            
            if i==1  % top row boundary
%                 u_Ne = 2*u_xstag(i-1,j) - u_xstag(i+1,j) ;
                u_Ne = 2*1 - 1*u_xstag(i,j) ; 

                u_ee = u_xstag(i,j+1) ;
                u_w = u_xstag(i,j-1) ;
                u_Se = u_xstag(i+1,j) ;
            elseif i==n_y  % bottom row boundary
%                 u_Se = 2*u_xstag(i,j) - u_xstag(i-1,j) ;
                u_Se = 2*0 - 1*u_xstag(i,j) ; 

                u_ee = u_xstag(i,j+1) ;
                u_Ne = u_xstag(i-1,j) ;
                u_w = u_xstag(i,j-1) ;
                           
            else            % inner rows
                u_ee = u_xstag(i,j+1) ;
                u_Ne = u_xstag(i-1,j) ;
                u_w = u_xstag(i,j-1) ;
                u_Se = u_xstag(i+1,j) ;
            end
 
            m_dot_e = rho*dy*w*u_xstag(i,j) ;
            m_dot_ee = rho*dy*w*u_xstag(i,j+1) ;
            m_dot_E = 0.5*(m_dot_e + m_dot_ee) ;
            
            m_dot_n = rho*dx*w*v_ystag(i,j-1) ;
            m_dot_nee = rho*dx*w*v_ystag(i,j) ;
            m_dot_ne = 0.5*(m_dot_n + m_dot_nee) ;

            m_dot_w = -rho*dy*w*u_xstag(i,j-1) ;
            m_dot_e = -rho*dy*w*u_xstag(i,j) ;
            m_dot_P = 0.5*(m_dot_w +m_dot_e) ;
            
            m_dot_s = -rho*dx*w*v_ystag(i+1,j-1) ;
            m_dot_see = -rho*dx*w*v_ystag(i+1,j) ;
            m_dot_se = 0.5*(m_dot_s + m_dot_see) ;

            
            a_ee_t = mu*dy*w/dx ;
            a_Ne_t = mu*dx*w/dy ;
            a_w_t = mu*dy*w/dx ;
            a_Se_t = mu*dx*w/dy ;

            a_ee_c = -m_dot_E/2 ;
            a_Ne_c = -m_dot_ne/2 ;
            a_w_c = -m_dot_P/2 ;
            a_Se_c = -m_dot_se/2 ;

            a_ee = a_ee_t + a_ee_c ;
            a_Ne = a_Ne_t + a_Ne_c ;
            a_w = a_w_t + a_w_c ;
            a_Se = a_Se_t + a_Se_c ;

            summ = u_ee*a_ee + u_Ne*a_Ne + u_w*a_w + u_Se*a_Se  ;

            u_xstag_pred(i,j) =  (P_E*A_E - P_P*A_P - summ)/ae ;

        end
    end



    
    % Starting with the Top-left-1 corner of yStag grid
    % we plug-in values at all points except the first&last rows of v_yStag matrix
    for i = 2:n_y               % inner row-wise             
        for j = 1:n_x            % column-wise

            % computing an coefficient at (i,j) coordinate in yStag grid
            an = get_an(i,j, u_xstag, v_ystag) ; 
            
            P_N = P_unstag(i-1,j) ;
            P_P = P_unstag(i,j) ;
            A_N = dx*w ;
            A_P = dx*w ;
            
            if j==1 % left column boundary
%                 v_nW = 2*v_ystag(i,j) - v_ystag(i,j+1) ;
                v_nW = 2*0 - 1*v_ystag(i,j) ;

                v_nE = v_ystag(i,j+1) ;
                v_Nn = v_ystag(i-1,j) ;
                v_s = v_ystag(i+1,j) ;
            elseif j==n_x % boundaries 
%                 v_nE = 2*v_ystag(i,j) - v_ystag(i,j-1) ;
                v_nE = 2*0 - 1*v_ystag(i,j) ;

                v_Nn = v_ystag(i-1,j) ;
                v_nW = v_ystag(i,j-1) ;
                v_s = v_ystag(i+1,j) ;
                           
            else
                v_nE = v_ystag(i,j+1) ;
                v_Nn = v_ystag(i-1,j) ;
                v_nW = v_ystag(i,j-1) ;
                v_s = v_ystag(i+1,j) ;
            end
 
            m_dot_Ne = rho*dy*w*u_xstag(i-1,j+1) ;
            m_dot_e = rho*dy*w*u_xstag(i,j+1) ;
            m_dot_ne = 0.5*(m_dot_Ne + m_dot_e) ;

            m_dot_n = rho*dx*w*v_ystag(i,j) ;
            m_dot_Nn = rho*dx*w*v_ystag(i-1,j) ;
            m_dot_N = 0.5*(m_dot_n + m_dot_Nn) ;
            
            m_dot_w = -rho*dy*w*u_xstag(i,j) ;
            m_dot_Nw = -rho*dy*w*u_xstag(i-1,j) ;
            m_dot_nw = 0.5*(m_dot_Nw +m_dot_w) ;
            
            m_dot_s = -rho*dx*w*v_ystag(i+1,j) ;
            m_dot_n = -rho*dx*w*v_ystag(i,j) ;
            m_dot_P = 0.5*(m_dot_n + m_dot_s) ;

            
            a_nE_t = mu*dy*w/dx ;
            a_Nn_t = mu*dx*w/dy ;
            a_nW_t = mu*dy*w/dx ;
            a_s_t = mu*dx*w/dy ;

            a_nE_c = -m_dot_ne/2 ;
            a_Nn_c = -m_dot_N/2 ;
            a_nW_c = -m_dot_nw/2 ;
            a_s_c = -m_dot_P/2 ;

            a_nE = a_nE_t + a_nE_c ;
            a_Nn = a_Nn_t + a_Nn_c ;
            a_nW = a_nW_t + a_nW_c ;
            a_s = a_s_t + a_s_c ;

            summ = v_nE*a_nE + v_Nn*a_Nn + v_nW*a_nW + v_s*a_s  ;

            v_ystag_pred(i,j) =  (P_N*A_N - P_P*A_P - summ)/an ;

        end
    end

    % Step3 : Find Pressure corrections   
    % according to the corrected velocity values
    % Starting with the Top-left corner of unstag grid
    P_unstag_c_old = P_unstag_c ;            % Just to ensure an explicit way of solving
    for i = 1:n_y         % row-wise             
        for j = 1:n_x     % column-wise

            u_e = u_xstag_pred(i,j+1) ;
            u_w = u_xstag_pred(i,j);
            A_e = dy*w ;
            A_w = dy*w ;

            v_n = v_ystag_pred(i,j) ;
            v_s = v_ystag_pred(i+1,j) ;
            A_n = dx*w ;
            A_s = dx*w ;

            b = A_e*u_e - A_w*u_w + A_n*v_n - A_s*v_s ;
            
            A_e = dy*w ;
            A_n = dx*w ;
            A_w = dy*w ;
            A_s = dx*w ;

            A_N = dx*w ;
            A_S = dx*w ;
            A_E = dy*w ;
            A_W = dy*w ;
            A_P = 0.5*(dx+dy)*w ;

            a_e = get_ae(i,j+1, u_xstag_pred, v_ystag_pred) ;
            a_w = get_ae(i,j, u_xstag_pred, v_ystag_pred) ;
            a_n = get_an(i,j, u_xstag_pred, v_ystag_pred) ;
            a_s = get_an(i+1,j, u_xstag_pred, v_ystag_pred) ;

            coeff_E = A_e*A_E/a_e ;
            coeff_N = A_n*A_N/a_n ;
            coeff_W = A_w*A_W/a_w ;
            coeff_S = A_s*A_S/a_s ;
           

            % finding P_c at neighbours of (i,j) in unstag grid
            if (i>1 && i<n_y) && (j>1 && j<n_x)
                % for inner grid points
                P_E_c = P_unstag_c_old(i,j+1) ;
                P_N_c = P_unstag_c_old(i-1,j) ;
                P_W_c = P_unstag_c_old(i,j-1) ;
                P_S_c = P_unstag_c_old(i+1,j) ;

                coeff_P = A_P*(A_e/a_e + A_w/a_w + A_n/a_n +A_s/a_s) ;
                P_unstag_c(i,j) = (coeff_E*P_E_c + coeff_N*P_N_c + coeff_W*P_W_c +coeff_S*P_S_c +b)/coeff_P ;
            elseif i == 1 && j==1 
                % for top-left corner
                P_E_c = P_unstag_c_old(i,j+1) ;
%                 P_N_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i+1,j) ;
%                 P_W_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i,j+1) ;
                P_S_c = P_unstag_c_old(i+1,j) ;
                
%               We try to leave out all the variables which came from the
%               vel corrections on the west and north sides
                coeff_P = A_P*(A_e/a_e +A_s/a_s) ;
                P_unstag_c(i,j) = (coeff_E*P_E_c +coeff_S*P_S_c +b)/coeff_P ;
            elseif i == 1 && (j > 1 && j<n_x)
                % for inner points on the top side
                P_E_c = P_unstag_c_old(i,j+1) ;
%                 P_N_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i+1,j) ;
                P_W_c = P_unstag_c_old(i,j-1) ;
                P_S_c = P_unstag_c_old(i+1,j) ;

                coeff_P = A_P*(A_e/a_e + A_w/a_w +A_s/a_s) ;
                P_unstag_c(i,j) = (coeff_E*P_E_c + coeff_W*P_W_c +coeff_S*P_S_c +b)/coeff_P ;
            elseif i ==1 && j== n_x 
                % for top-right corner
                P_E_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i,j-1) ;
                P_N_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i+1,j) ;
                P_W_c = P_unstag_c_old(i,j-1) ;
                P_S_c = P_unstag_c_old(i+1,j) ;

                coeff_P = A_P*(A_w/a_w +A_s/a_s) ;
                P_unstag_c(i,j) = ( coeff_W*P_W_c +coeff_S*P_S_c +b)/coeff_P ;
            elseif (i > 1 && i < n_y) && j ==n_x
                % for inner points on the right side
                P_E_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i,j-1) ;
                P_N_c = P_unstag_c_old(i-1,j) ;
                P_W_c = P_unstag_c_old(i,j-1) ;
                P_S_c = P_unstag_c_old(i+1,j) ;

                coeff_P = A_P*( A_w/a_w + A_n/a_n +A_s/a_s) ;
                P_unstag_c(i,j) = ( coeff_N*P_N_c + coeff_W*P_W_c +coeff_S*P_S_c +b)/coeff_P ;
            elseif i==n_y && j==n_x
                % for bottom-right corner
                P_E_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i,j-1) ;
                P_N_c = P_unstag_c_old(i-1,j) ;
                P_W_c = P_unstag_c_old(i,j-1) ;
                P_S_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i-1,j) ;  

                coeff_P = A_P*( A_w/a_w + A_n/a_n ) ;
                P_unstag_c(i,j) = ( coeff_N*P_N_c + coeff_W*P_W_c +b)/coeff_P ;   
            elseif i==n_y && (j>1 && j< n_x)
                % for inner points on the bottom side
                P_E_c = P_unstag_c_old(i,j+1) ;
                P_N_c = P_unstag_c_old(i-1,j) ;
                P_W_c = P_unstag_c_old(i,j-1) ;
                P_S_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i-1,j) ;    

                coeff_P = A_P*(A_e/a_e + A_w/a_w + A_n/a_n ) ;
                P_unstag_c(i,j) = (coeff_E*P_E_c + coeff_N*P_N_c + coeff_W*P_W_c  +b)/coeff_P ;

% Note : We leave out pressure correction at the bottom left pressure point
% Through this we ensure that P_at bottom left is equal to 1
            elseif i==n_y && j==1
                % for bottom-left corner
                P_E_c = P_unstag_c_old(i,j+1) ;
                P_N_c = P_unstag_c_old(i-1,j) ;
                P_W_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i,j+1) ;
                P_S_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i-1,j) ; 

                coeff_P = A_P*(A_e/a_e + A_n/a_n ) ;
                P_unstag_c(i,j) = (coeff_E*P_E_c + coeff_N*P_N_c +b)/coeff_P ;

            elseif (i>1 && i<n_y) && j==1
                % for inner points on the left side
                P_E_c = P_unstag_c_old(i,j+1) ;
                P_N_c = P_unstag_c_old(i-1,j) ;
                P_W_c = 2*P_unstag_c_old(i,j) - P_unstag_c_old(i,j+1) ;
                P_S_c = P_unstag_c_old(i+1,j) ; 

                coeff_P = A_P*(A_e/a_e + A_n/a_n +A_s/a_s) ;
                P_unstag_c(i,j) = (coeff_E*P_E_c + coeff_N*P_N_c +coeff_S*P_S_c +b)/coeff_P ;
            end
        


        end
    end


    % Step4 : Find Velocity corrections

    % Starting with the Top-left-1 corner of xStag grid
    % we leave out the first & last columns since u is known there
    for i = 1:n_y               % row-wise             
        for j = 2:n_x               % inner column-wise

            ae = get_ae(i,j,u_xstag_pred, v_ystag_pred) ;

            A_E = dy*w ;
            A_P = dy*w ;

            P_E_c = P_unstag_c(i,j) ;
            P_P_c = P_unstag_c(i,j-1) ;

            u_xstag_c(i,j) = (P_E_c*A_E - P_P_c*A_P)/ae ;
        end
    end

    % Starting with the Top-left-1 corner of yStag grid
    % we leave out the first & last rows since v is know there
    for i = 2:n_y           % inner row-wise             
        for j = 1:n_x               % column-wise
            
            an = get_an(i,j, u_xstag_pred, v_ystag_pred) ;

            P_N_c = P_unstag_c(i-1,j) ;
            P_P_c = P_unstag_c(i,j);

            A_P = dx*w ;
            A_N = dx*w ;

            v_ystag_c(i,j) = (P_N_c*A_N - P_P_c*A_P)/an ;
        end
    end


    % Step5 : Correct the Pressure and Velocities

    % Correct pressure starting from the top left corner
    for i = 1:n_y
        for j = 1:n_x 
            P_unstag_new(i,j) = P_unstag(i,j) + alpha_P*P_unstag_c(i,j) ;
        end
    end

    % Correct u starting from the top left corner
    % we leave out the first & last columns since u is fixed there
    for i = 1:n_y
        for j = 2:n_x 
            u_xstag_new(i,j) = u_xstag_pred(i,j) + alpha_v*u_xstag_c(i,j) ;
        end
    end

    % Correct v starting from the top left corner
    % we leave out the first & last rows since v is fixed there
    for i = 2:n_y
        for j = 1:n_x 
            v_ystag_new(i,j) = v_ystag_pred(i,j) + alpha_v*v_ystag_c(i,j) ;
        end
    end


    % Step6 : Check for Convergence

    % We find the LHS of the continuity eqn satisfied by the (u,v,P) at
    % every internal grid point of unstag grid
    LHS = ones(n_y,n_x) ;

    % Starting with the Top-left corner of unstag grid
    for i = 1:n_y           % row-wise             
        for j = 1:n_x       % column-wise

            u_e = u_xstag_new(i,j+1) ;
            u_w = u_xstag_new(i,j);
            A_e = dy*w ;
            A_w = dy*w ;

            v_n = v_ystag_new(i,j) ;
            v_s = v_ystag_new(i+1,j) ;
            A_n = dx*w ;
            A_s = dx*w ;

            LHS(i,j) = A_e*u_e - A_w*u_w + A_n*v_n - A_s*v_s ;
        end
    end
    
    error_val = max(max(abs(LHS))) ;
    str_disp = sprintf("iter: %f error %f", iter, error_val) ;
    disp(str_disp) ;

    if error_val <= eps 
        disp("---------------------------------------------------------") ;
        str_disp = sprintf("The solution Converged at iter: %f", iter) ;
        disp(str_disp) ;
        break
    end
    
    % re-initialise 
    u_xstag = u_xstag_new ;
    v_ystag = v_ystag_new ;
    P_unstag = P_unstag_new ;


end



%% Plotting the Contour and Velocity at some lines
% Beautify Contour and Velocity Plots

% Define a colormap for all plots
% Define a colormap for all plots
colormap_turbo = turbo; % High-contrast, perceptually uniform colormap

% Plot the u-velocity contour
figure(1);
u_xstag_averaged = zeros(n_y, n_x);
for j = 1:n_x
    u_xstag_averaged(:,j) = 0.5 * (u_xstag_new(:, j) + u_xstag_new(:, j + 1));
end
contourf(flip(u_xstag_averaged, 1), 20, 'LineStyle', 'none'); % Use 20 contour levels
colormap(colormap_turbo);
colorbar; 
caxis([-max(abs(u_xstag_averaged(:))), max(abs(u_xstag_averaged(:)))]); % Symmetric color scale
xlabel('X-direction', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y-direction', 'FontSize', 12, 'FontWeight', 'bold');
title('u-Velocity Contour Plot', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot the v-velocity contour
figure(2);
v_ystag_averaged = zeros(n_y, n_x);
for i = 1:n_y
    v_ystag_averaged(i,:) = 0.5 * (v_ystag_new(i,:) + v_ystag_new(i + 1,:));
end
contourf(flip(v_ystag_averaged, 1), 20, 'LineStyle', 'none'); 
colormap(colormap_turbo);
colorbar;
caxis([-max(abs(v_ystag_averaged(:))), max(abs(v_ystag_averaged(:)))]); 
xlabel('X-direction', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y-direction', 'FontSize', 12, 'FontWeight', 'bold');
title('v-Velocity Contour Plot', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot the Pressure contour
figure(3);
P_unstag_reduced = P_unstag_new(:, :);
contourf(flip(P_unstag_reduced, 1), 20, 'LineStyle', 'none'); 
colormap(colormap_turbo);
colorbar;
xlabel('X-direction', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y-direction', 'FontSize', 12, 'FontWeight', 'bold');
title('Pressure Contour Plot', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot the u-velocity along x = 0.5 line
figure(4);
u_along_x = flip(u_xstag_averaged(:, floor(n_x/2)), 1);
plot(linspace(0,1,n_y), u_along_x, 'LineWidth', 1.5, 'Color', 'b');
xlabel('Position along x = 0.5 line', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('u-Velocity', 'FontSize', 12, 'FontWeight', 'bold');
title('u-Velocity Along x = 0.5', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot the v-velocity along y = 0.5 line
figure(5);
v_along_y = v_ystag_averaged(floor(n_y/2), :);
plot(linspace(0,1,n_x), v_along_y, 'LineWidth', 1.5, 'Color', 'r');
xlabel('Position along y = 0.5 line', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('v-Velocity', 'FontSize', 12, 'FontWeight', 'bold');
title('v-Velocity Along y = 0.5', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot the Pressure along x = 0.5 line
figure(6);
P_along_x = flip(P_unstag_new(:, floor(n_x/2)), 1);
plot(linspace(0,1,n_y), P_along_x, 'LineWidth', 1.5, 'Color', 'g');
xlabel('Position along x = 0.5 line', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Pressure', 'FontSize', 12, 'FontWeight', 'bold');
title('Pressure Along x = 0.5', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot the Pressure along y = 0.5 line
figure(7);
P_along_y = P_unstag_new(floor(n_y/2), :);
plot(linspace(0,1,n_x), P_along_y, 'LineWidth', 1.5, 'Color', 'm');
xlabel('Position along y = 0.5 line', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Pressure', 'FontSize', 12, 'FontWeight', 'bold');
title('Pressure Along y = 0.5', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot the velocity streamlines
figure(8);
flipped_u = flip(u_xstag_averaged, 1);
flipped_v = flip(v_ystag_averaged, 1);
[X, Y] = meshgrid(linspace(0, L_x, n_x), linspace(0, L_y, n_y));
streamslice(X, Y, flipped_u, flipped_v, 'ArrowStyle', 'classic'); % Classic arrow style
xlabel('X-direction', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y-direction', 'FontSize', 12, 'FontWeight', 'bold');
title('Velocity Streamlines', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
axis tight;





%% Useful Coefficient Formula's

% return coeff ae for e points at (i,j) in the xStag grid
function ae = get_ae(i,j, u_xstag, v_ystag)
    % Note: m_dot at up and down sides here will take the u velocity 
    % Note: m_dot at left and right sides here will take the v velocity (from v_ystag)


    global n_x
    global n_y
    global dx
    global dy
    global w
    global rho
    global mu 

    if (i > 1 && i < n_y) && (j > 1 && j < (n_x+1)) 
        % for inner grid points
        m_dot_e = rho*dy*w*u_xstag(i,j) ;
        m_dot_ee = rho*dy*w*u_xstag(i,j+1) ;
        m_dot_E = 0.5*(m_dot_e + m_dot_ee) ;
        
        m_dot_n = rho*dx*w*v_ystag(i,j-1) ;
        m_dot_nee = rho*dx*w*v_ystag(i,j) ;
        m_dot_ne = 0.5*(m_dot_n + m_dot_nee);

        m_dot_e = -rho*dy*w*u_xstag(i,j) ; % is not the same as before
        m_dot_w = -rho*dy*w*u_xstag(i,j-1) ;
        m_dot_P = 0.5*(m_dot_w + m_dot_e) ;

        m_dot_s = -rho*dx*w*v_ystag(i+1,j-1) ;
        m_dot_seE = -rho*dx*w*v_ystag(i+1,j) ;
        m_dot_se = 0.5*(m_dot_s + m_dot_seE) ;

        ae_c = -0.5*(m_dot_E + m_dot_ne + m_dot_P + m_dot_se) ;

    elseif (i == 1 && j ==1 )
        % For top-left corner
        
        m_dot_e = rho*dy*w*u_xstag(i,j) ; % equals 0 since we have u=0 permenant on the left side
        m_dot_ee = rho*dy*w*u_xstag(i,j+1) ;
        m_dot_E = 0.5*(m_dot_e + m_dot_ee) ;
        
%         m_dot_nee = rho*dx*w*v_ystag(i,j) ;
%         m_dot_neee = rho*dx*w*v_ystag(i,j+1) ;
%         m_dot_ne = 1.5*m_dot_nee - 0.5*m_dot_neee ;
        m_dot_ne = rho*dx*w*0 ; % manually setting it to 0 since we are not ensuring v_left = 0 anywhere else

        m_dot_e = -rho*dy*w*u_xstag(i,j) ; % is not the same as before
        m_dot_ee = -rho*dy*w*u_xstag(i,j+1) ;
        m_dot_P = 1.5*m_dot_e - 0.5*m_dot_ee ;

%         m_dot_seE = -rho*dx*w*v_ystag(i+1,j) ;
%         m_dot_seEE = -rho*dx*w*v_ystag(i+1,j+1) ;
%         m_dot_se = 1.5*m_dot_seE - 0.5*m_dot_seEE ;
        m_dot_se = rho*dx*w*0 ; % manually setting it to zero

        ae_c = -0.5*(m_dot_E + m_dot_ne + m_dot_P + m_dot_se) ;

    elseif (i == 1 ) && ( j > 1 && j < (n_x+1) ) 
        % For inner part of top-side
        m_dot_e = rho*dy*w*u_xstag(i,j) ; 
        m_dot_ee = rho*dy*w*u_xstag(i,j+1) ;
        m_dot_E = 0.5*(m_dot_e + m_dot_ee) ;
        
%         m_dot_n = rho*dx*w*v_ystag(i,j-1) ;
%         m_dot_nee = rho*dx*w*v_ystag(i,j) ;
%         m_dot_ne = 1.5*m_dot_n - 0.5*m_dot_nww ;
        m_dot_ne = rho*dx*w*0 ;   % manually setting it to zero

        m_dot_e = -rho*dy*w*u_xstag(i,j) ; 
        m_dot_w = -rho*dy*w*u_xstag(i,j-1) ;
        m_dot_P = 0.5*(m_dot_w + m_dot_e) ;

        m_dot_s = -rho*dx*w*v_ystag(i+1,j-1) ;
        m_dot_seE = -rho*dx*w*v_ystag(i+1,j) ;
        m_dot_se = 0.5*(m_dot_s + m_dot_seE) ;

        ae_c = -0.5*(m_dot_E + m_dot_ne + m_dot_P + m_dot_se) ;


    elseif (i == 1 && i < n_y) && j ==(n_x+1)
        % For top-right corner
        m_dot_e = rho*dy*w*u_xstag(i,j) ;   % equals 0
        m_dot_w = rho*dy*w*u_xstag(i,j-1) ;
        m_dot_E = 1.5*m_dot_e - 0.5*m_dot_w ;
        
%         m_dot_n = rho*dx*w*v_ystag(i,j-1) ;
%         m_dot_nww = rho*dx*w*v_ystag(i,j-2) ;
%         m_dot_ne = 1.5*m_dot_n - 0.5*m_dot_nww ;
        m_dot_ne = rho*dx*w*0 ; % manually setting it to zero

        m_dot_e = -rho*dy*w*u_xstag(i,j) ; % is not the same as before
        m_dot_w = -rho*dy*w*u_xstag(i,j-1) ;
        m_dot_P = 0.5*(m_dot_w + m_dot_e) ;

%         m_dot_s = -rho*dx*w*v_ystag(i+1,j-1) ;
%         m_dot_sWw = -rho*dx*w*v_ystag(i+1,j-2) ;
%         m_dot_se = 1.5*m_dot_s - 0.5*m_dot_sWw ;
        m_dot_se = rho*dx*w*0 ; % manually setting it to zero

        ae_c = -0.5*(m_dot_E + m_dot_ne + m_dot_P + m_dot_se) ;

    elseif (i > 1 && i < n_y) && j ==(n_x+1)
        % For inner part of right side
        m_dot_e = rho*dy*w*u_xstag(i,j) ;  
        m_dot_w = rho*dy*w*u_xstag(i,j-1) ;
        m_dot_E = 1.5*m_dot_e - 0.5*m_dot_w ;
        
%         m_dot_n = rho*dx*w*v_ystag(i,j-1) ;
%         m_dot_nww = rho*dx*w*v_ystag(i,j-2) ;
%         m_dot_ne = 1.5*m_dot_n - 0.5*m_dot_nww ;
        m_dot_ne = rho*dx*w*0 ; % manually setting it to zero

        m_dot_e = -rho*dy*w*u_xstag(i,j) ; % is not the same as before
        m_dot_w = -rho*dy*w*u_xstag(i,j-1) ;
        m_dot_P = 0.5*(m_dot_w + m_dot_e) ;

%         m_dot_s = -rho*dx*w*v_ystag(i+1,j-1) ;
%         m_dot_sWw = -rho*dx*w*v_ystag(i+1,j-2) ;
%         m_dot_se = 1.5*m_dot_s - 0.5*m_dot_sWw ;
        m_dot_se = rho*dx*w*0 ; % manually setting it to zero

        ae_c = -0.5*(m_dot_E + m_dot_ne + m_dot_P + m_dot_se) ;

    elseif i== n_y && j ==(n_x+1)
        % For bottom right corner
        m_dot_e = rho*dy*w*u_xstag(i,j) ;  
        m_dot_w = rho*dy*w*u_xstag(i,j-1) ;
        m_dot_E = 1.5*m_dot_e - 0.5*m_dot_w ;
        
%         m_dot_n = rho*dx*w*v_ystag(i,j-1) ;
%         m_dot_nww = rho*dx*w*v_ystag(i,j-2) ;
%         m_dot_ne = 1.5*m_dot_n - 0.5*m_dot_nww ;
        m_dot_ne = rho*dx*w*0 ; % manually setting it to zero

        m_dot_e = -rho*dy*w*u_xstag(i,j) ; % is not the same as before
        m_dot_w = -rho*dy*w*u_xstag(i,j-1) ;
        m_dot_P = 0.5*(m_dot_w + m_dot_e) ;

%         m_dot_s = -rho*dx*w*v_ystag(i+1,j-1) ;
%         m_dot_sWw = -rho*dx*w*v_ystag(i+1,j-2) ;
%         m_dot_se = 1.5*m_dot_s - 0.5*m_dot_sWw ;
        m_dot_se = rho*dx*w*0 ; % manually setting it to zero

        ae_c = -0.5*(m_dot_E + m_dot_ne + m_dot_P + m_dot_se) ;

    elseif i== n_y && ( j>1 && j < (n_x+1) )
        % For inner part of bottom side
        m_dot_e = rho*dy*w*u_xstag(i,j) ;  
        m_dot_ee = rho*dy*w*u_xstag(i,j+1) ;
        m_dot_E = 0.5*(m_dot_e + m_dot_ee) ;
        
        m_dot_n = rho*dx*w*v_ystag(i,j-1) ;
        m_dot_nee = rho*dx*w*v_ystag(i,j) ;
        m_dot_ne = 0.5*(m_dot_n + m_dot_nee) ;

        m_dot_e = -rho*dy*w*u_xstag(i,j) ; % is not the same as before
        m_dot_w = -rho*dy*w*u_xstag(i,j-1) ;
        m_dot_P = 0.5*(m_dot_w + m_dot_e) ;

        m_dot_s = -rho*dx*w*v_ystag(i+1,j-1) ;   % equals zero
        m_dot_seE = -rho*dx*w*v_ystag(i+1,j) ;   % equals zero
        m_dot_se = 0.5*(m_dot_s + m_dot_seE) ; 

        ae_c = -0.5*(m_dot_E + m_dot_ne + m_dot_P + m_dot_se) ;
    
    elseif i== n_y && j==1

        % For bottom left corner
        m_dot_e = rho*dy*w*u_xstag(i,j) ;  
        m_dot_ee = rho*dy*w*u_xstag(i,j+1) ;
        m_dot_E = 0.5*(m_dot_e + m_dot_ee) ;
        
%         m_dot_nee = rho*dx*w*v_ystag(i,j) ;
%         m_dot_neee = rho*dx*w*v_ystag(i,j+1) ;
%         m_dot_ne = 1.5*m_dot_nee - 0.5*m_dot_neee ;
        m_dot_ne = rho*dx*w*0 ; % manually setting it to zero

        m_dot_e = -rho*dy*w*u_xstag(i,j) ; % is not the same as before
        m_dot_ee = -rho*dy*w*u_xstag(i,j+1) ;
        m_dot_P = 1.5*m_dot_e - 0.5*m_dot_ee ;


%         m_dot_seE = -rho*dx*w*v_ystag(i+1,j) ;
%         m_dot_seEE = -rho*dx*w*v_ystag(i+1,j+1) ;
%         m_dot_se = 1.5*m_dot_seE - 0.5*m_dot_seEE ;
        m_dot_se = rho*dx*w*0 ; % manually setting it to zero

        ae_c = -0.5*(m_dot_E + m_dot_ne + m_dot_P + m_dot_se) ;

    elseif (i>1 && i<n_y) && j==1

        % For inner part of left side
        m_dot_e = rho*dy*w*u_xstag(i,j) ;  
        m_dot_ee = rho*dy*w*u_xstag(i,j+1) ;
        m_dot_E = 0.5*(m_dot_e + m_dot_ee) ;
        
%         m_dot_nee = rho*dx*w*v_ystag(i,j) ;
%         m_dot_neee = rho*dx*w*v_ystag(i,j+1) ;
%         m_dot_ne = 1.5*m_dot_nee - 0.5*m_dot_neee ;
        m_dot_ne = rho*dx*w*0 ; % manually setting it to zero

        m_dot_e = -rho*dy*w*u_xstag(i,j) ; % is not the same as before
        m_dot_ee = -rho*dy*w*u_xstag(i,j+1) ;
        m_dot_P = 1.5*m_dot_e - 0.5*m_dot_ee ;

%         m_dot_seE = -rho*dx*w*v_ystag(i+1,j) ;
%         m_dot_seEE = -rho*dx*w*v_ystag(i+1,j+1) ;
%         m_dot_se = 1.5*m_dot_seE - 0.5*m_dot_seEE ;
        m_dot_se = rho*dx*w*0 ; % manually setting it to zero

        ae_c = -0.5*(m_dot_E + m_dot_ne + m_dot_P + m_dot_se) ;

    end

    ae_t = -2*mu*( (dy*w)/dx + (dx*w)/dy ) ;

    ae =  ae_c + ae_t ;

end

























% return coeff 'an' for point at (i,j) in the yStag grid
function an = get_an(i,j, u_xstag, v_ystag )
    % Note: m_dot at up and down sides here will take the u velocity 
    % Note: m_dot at left and right sides here will take the v velocity (from v_ystag)
    
    global n_x
    global n_y
    global dx
    global dy
    global w
    global rho
    global mu 

    if ( i > 1 && i < (n_y+1) ) && (j > 1 && j < n_x) 
        % for inner grid points
        m_dot_e = rho*dy*w*u_xstag(i,j+1) ;
        m_dot_Ne = rho*dy*w*u_xstag(i-1,j+1) ;
        m_dot_ne = 0.5*(m_dot_e + m_dot_Ne) ;
        
        m_dot_n = rho*dx*w*v_ystag(i,j) ;
        m_dot_nn = rho*dx*w*v_ystag(i-1,j) ;
        m_dot_N = 0.5*(m_dot_n + m_dot_nn);

        m_dot_Nw = -rho*dy*w*u_xstag(i-1,j) ; 
        m_dot_w = -rho*dy*w*u_xstag(i,j) ;
        m_dot_nw = 0.5*(m_dot_w + m_dot_Nw) ;

        m_dot_n = -rho*dx*w*v_ystag(i,j) ;   % is not the same as before
        m_dot_s = -rho*dx*w*v_ystag(i+1,j) ;
        m_dot_P = 0.5*(m_dot_s + m_dot_n) ;

        an_c = -0.5*(m_dot_ne + m_dot_N  + m_dot_nw + m_dot_P) ;

    elseif (i == 1 && j ==1 )
        % For top-left corner
        
%         m_dot_e = rho*dy*w*u_xstag(i,j+1) ;
%         m_dot_Se = rho*dy*w*u_xstag(i+1,j+1) ;
%         m_dot_ne = 1.5*m_dot_e - 0.5*m_dot_Se ;
        m_dot_ne = rho*dy*w*1 ;   % manually setting u = 1 since we are not putting the u_top=1 boundary anywhere else

        m_dot_n = rho*dx*w*v_ystag(i,j) ;
        m_dot_s = rho*dx*w*v_ystag(i+1,j) ;
        m_dot_N = 1.5*m_dot_n - 0.5*m_dot_s;

%         m_dot_Sw = -rho*dy*w*u_xstag(i+1,j) ; 
%         m_dot_w = -rho*dy*w*u_xstag(i,j) ;
%         m_dot_nw = 1.5*m_dot_w - 0.5*m_dot_Sw ;
        m_dot_nw = -rho*dy*w*1 ;   % manually setting u = 1

        m_dot_n = -rho*dx*w*v_ystag(i,j) ;   % is not the same as before
        m_dot_s = -rho*dx*w*v_ystag(i+1,j) ;
        m_dot_P = 0.5*(m_dot_s + m_dot_n) ;

        an_c = -0.5*(m_dot_ne + m_dot_N  + m_dot_nw + m_dot_P) ;

    elseif (i == 1 ) && (j > 1 && j < n_x) 
        % For inner part of top-side
%         m_dot_e = rho*dy*w*u_xstag(i,j+1) ;
%         m_dot_Se = rho*dy*w*u_xstag(i+1,j+1) ;
%         m_dot_ne = 1.5*m_dot_e - 0.5*m_dot_Se ;
        m_dot_ne = rho*dy*w*1 ;   % manually setting u = 1
        
        m_dot_n = rho*dx*w*v_ystag(i,j) ;
        m_dot_s = rho*dx*w*v_ystag(i+1,j) ;
        m_dot_N = 1.5*m_dot_n - 0.5*m_dot_s;

%         m_dot_Sw = -rho*dy*w*u_xstag(i+1,j) ; 
%         m_dot_w = -rho*dy*w*u_xstag(i,j) ;
%         m_dot_nw = 1.5*m_dot_w - 0.5*m_dot_Sw ;
        m_dot_nw = -rho*dy*w*1 ;   % manually setting u = 1

        m_dot_n = -rho*dx*w*v_ystag(i,j) ;   % is not the same as before
        m_dot_s = -rho*dx*w*v_ystag(i+1,j) ;
        m_dot_P = 0.5*(m_dot_s + m_dot_n) ;

        an_c = -0.5*(m_dot_ne + m_dot_N  + m_dot_nw + m_dot_P) ;


    elseif ( i == 1 && j ==n_x )
        % For top-right corner
%         m_dot_e = rho*dy*w*u_xstag(i,j+1) ;
%         m_dot_Se = rho*dy*w*u_xstag(i+1,j+1) ;
%         m_dot_ne = 1.5*m_dot_e - 0.5*m_dot_Se ;
        m_dot_ne = rho*dy*w*1 ;   % manually setting u = 1
        
        m_dot_n = rho*dx*w*v_ystag(i,j) ;
        m_dot_s = rho*dx*w*v_ystag(i+1,j) ;
        m_dot_N = 1.5*m_dot_n - 0.5*m_dot_s;

%         m_dot_Sw = -rho*dy*w*u_xstag(i+1,j) ; 
%         m_dot_w = -rho*dy*w*u_xstag(i,j) ;
%         m_dot_nw = 1.5*m_dot_w - 0.5*m_dot_Sw ;
        m_dot_nw = -rho*dy*w*1 ;   % manually setting u = 1

        m_dot_n = -rho*dx*w*v_ystag(i,j) ;   % is not the same as before
        m_dot_s = -rho*dx*w*v_ystag(i+1,j) ;
        m_dot_P = 0.5*(m_dot_s + m_dot_n) ;

        an_c = -0.5*(m_dot_ne + m_dot_N  + m_dot_nw + m_dot_P) ;

    elseif ( i > 1 && i < (n_y+1) ) && j == n_x
        % For inner part of right side
        m_dot_e = rho*dy*w*u_xstag(i,j+1) ;  % equals 0
        m_dot_Ne = rho*dy*w*u_xstag(i-1,j+1) ; % equals 0
        m_dot_ne = 0.5*(m_dot_e + m_dot_Ne) ;
        
        m_dot_n = rho*dx*w*v_ystag(i,j) ;
        m_dot_nn = rho*dx*w*v_ystag(i-1,j) ;
        m_dot_N = 0.5*(m_dot_n + m_dot_nn);

        m_dot_Nw = -rho*dy*w*u_xstag(i-1,j) ; 
        m_dot_w = -rho*dy*w*u_xstag(i,j) ;
        m_dot_nw = 0.5*(m_dot_w + m_dot_Nw) ;

        m_dot_n = -rho*dx*w*v_ystag(i,j) ;   % is not the same as before
        m_dot_s = -rho*dx*w*v_ystag(i+1,j) ;
        m_dot_P = 0.5*(m_dot_s + m_dot_n) ;

        an_c = -0.5*(m_dot_ne + m_dot_N  + m_dot_nw + m_dot_P) ;

    elseif (i == (n_y+1) && j == n_x )
        % For bottom right corner
%         m_dot_NNe = rho*dy*w*u_xstag(i-2,j+1) ;
%         m_dot_Ne = rho*dy*w*u_xstag(i-1,j+1) ;
%         m_dot_ne = 1.5*m_dot_Ne - 0.5*m_dot_NNe ;
        m_dot_ne = rho*dy*w*0 ;   % manually setting u = 0
        
        m_dot_n = rho*dx*w*v_ystag(i,j) ;
        m_dot_nn = rho*dx*w*v_ystag(i-1,j) ;
        m_dot_N = 0.5*(m_dot_n + m_dot_nn);

%         m_dot_Nw = -rho*dy*w*u_xstag(i-1,j) ; 
%         m_dot_NNw = -rho*dy*w*u_xstag(i-2,j) ;
%         m_dot_nw = 1.5*m_dot_Nw - 0.5*m_dot_NNw ;
        m_dot_nw = -rho*dy*w*0 ;   % manually setting u = 0

        m_dot_n = -rho*dx*w*v_ystag(i,j) ;   % is not the same as before
        m_dot_Nn = -rho*dx*w*v_ystag(i-1,j) ;
        m_dot_P = 1.5*m_dot_n - 0.5*m_dot_Nn ;

        an_c = -0.5*(m_dot_ne + m_dot_N  + m_dot_nw + m_dot_P) ;

    elseif i== (n_y+1) && ( j>1 && j <n_x )
        % For inner part of bottom side
%         m_dot_NNe = rho*dy*w*u_xstag(i-2,j+1) ;
%         m_dot_Ne = rho*dy*w*u_xstag(i-1,j+1) ;
%         m_dot_ne = 1.5*m_dot_Ne - 0.5*m_dot_NNe ;
        m_dot_ne = rho*dy*w*0 ;   % manually setting u = 0
        
        m_dot_n = rho*dx*w*v_ystag(i,j) ;
        m_dot_nn = rho*dx*w*v_ystag(i-1,j) ;
        m_dot_N = 0.5*(m_dot_n + m_dot_nn);

%         m_dot_Nw = -rho*dy*w*u_xstag(i-1,j) ; 
%         m_dot_NNw = -rho*dy*w*u_xstag(i-2,j) ;
%         m_dot_nw = 1.5*m_dot_Nw - 0.5*m_dot_NNw ;
        m_dot_nw = -rho*dy*w*0 ;   % manually setting u = 0

        m_dot_n = -rho*dx*w*v_ystag(i,j) ;   % is not the same as before
        m_dot_Nn = -rho*dx*w*v_ystag(i-1,j) ;
        m_dot_P = 1.5*m_dot_n - 0.5*m_dot_Nn ;

        an_c = -0.5*(m_dot_ne + m_dot_N  + m_dot_nw + m_dot_P) ;
    
    elseif i== (n_y+1) && j==1

        % For bottom left corner
%         m_dot_NNe = rho*dy*w*u_xstag(i-2,j+1) ;
%         m_dot_Ne = rho*dy*w*u_xstag(i-1,j+1) ;
%         m_dot_ne = 1.5*m_dot_Ne - 0.5*m_dot_NNe ;
        m_dot_ne = rho*dy*w*0 ;   % manually setting u = 0
        
        m_dot_n = rho*dx*w*v_ystag(i,j) ;
        m_dot_nn = rho*dx*w*v_ystag(i-1,j) ;
        m_dot_N = 0.5*(m_dot_n + m_dot_nn);

%         m_dot_Nw = -rho*dy*w*u_xstag(i-1,j) ; 
%         m_dot_NNw = -rho*dy*w*u_xstag(i-2,j) ;
%         m_dot_nw = 1.5*m_dot_Nw - 0.5*m_dot_NNw ;
        m_dot_nw = -rho*dy*w*0 ;   % manually setting u = 0

        m_dot_n = -rho*dx*w*v_ystag(i,j) ;   % is not the same as before
        m_dot_Nn = -rho*dx*w*v_ystag(i-1,j) ;
        m_dot_P = 1.5*m_dot_n - 0.5*m_dot_Nn ;

        an_c = -0.5*(m_dot_ne + m_dot_N  + m_dot_nw + m_dot_P) ;

    elseif ( i>1 && i<(n_y+1) ) && j==1

        % For inner part of left side
        m_dot_e = rho*dy*w*u_xstag(i,j+1) ;
        m_dot_Ne = rho*dy*w*u_xstag(i-1,j+1) ;
        m_dot_ne = 0.5*(m_dot_e + m_dot_Ne) ;
        
        m_dot_n = rho*dx*w*v_ystag(i,j) ;
        m_dot_nn = rho*dx*w*v_ystag(i-1,j) ;
        m_dot_N = 0.5*(m_dot_n + m_dot_nn);

        m_dot_Nw = -rho*dy*w*u_xstag(i-1,j) ; 
        m_dot_w = -rho*dy*w*u_xstag(i,j) ;
        m_dot_nw = 0.5*(m_dot_w + m_dot_Nw) ;

        m_dot_n = -rho*dx*w*v_ystag(i,j) ;   % is not the same as before
        m_dot_s = -rho*dx*w*v_ystag(i+1,j) ;
        m_dot_P = 0.5*(m_dot_s + m_dot_n) ;

        an_c = -0.5*(m_dot_ne + m_dot_N  + m_dot_nw + m_dot_P) ;

    end


    an_t = -2*mu*( (dy*w)/dx + (dx*w)/dy ) ;

    an = an_c + an_t ;

end














