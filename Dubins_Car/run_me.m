% add path to Casadi
% addpath('C:\Users\wiles\Downloads\casadi-3.6.7-windows64-matlab2018b');
compute_jacobians = 1;

if(compute_jacobians)
    clearvars -except compute_jacobians
else
    clearvars -except Jacobian_u_G Jacobian_u_J Rho_a_sym Jacobian_Rho_a compute_jacobians Rho_a_final_step Jacobian_Rho_a_final_step
end

import casadi.*

verbose = 0; % 1 for more info

obstacle_radius = 1;
Q1 = (1/obstacle_radius^2)*eye(2);

Q = [Q1, Q1];
c = [[1;2],[4-0.95;2]];

N = 20; % horizon length
color = 'ro-';
max_steps = 50;
ds = 0.25; % step size in predictor-corrector
min_car_turning_radius = 1/6;
v = 1; % car's constant speed
x1 = [0;1;-pi/2];
x_target = [2;3;0];

% many_obs_flag = 1;
% 
% if(many_obs_flag)
%     Q = [];
%     for j = 1:10
%         Q = [Q, Q1];
%     end
%     c = [[1;2.5],[2;0], [4;2], [5;-3], [6;-1], [7;1], [8;-1], [9;1], [10;3], [11;0]];
% 
%     N = 46; % horizon length
%     color = 'ro-';
%     max_steps = 150;
%     ds = 0.25;
%     min_car_turning_radius = 1/8;
%     v = 3;
%     x1 = [0;0;0];
%     x_target = [13;2;0];
% else
% 
% end

figure(1);
plot_ellipsoids(Q, c)
plot(x1(1), x1(2), 'square')
plot(x_target(1), x_target(2), 'kx')
xlabel('$x_1$', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);
ylabel('$x_2$', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);

% control's Euclidean norm is constrained:
Q_c = eye(1);

n = 3; % state dimension
m = 1; % control dimension
p = size(Q,2)/2; % number of state constraints
q = 1; % number of control constraints

r = N*m;
s = N*(p + q);

lambda = SX.sym('lambda',1);
u = SX.sym('u', r);
mu = SX.sym('mu', s);

% select a random a
u0 = rand(r,1);
b0 = get_b0(u0,x1,N,m,r,s,Q,c,Q_c,min_car_turning_radius,v);

if(G(0,u0,x1,N,m,r,Q,c,Q_c,min_car_turning_radius,v) < b0)
    if(verbose)
        disp('b0 is good!');
    end
end

c0 = ones(s,1);         % s
a = [u0;b0;c0]; % should be (r + 2*s) x 1 

% plot initial guess
x_tot = [];
for k=1:N+1
    x_tot = [x_tot, sol_diff_eqn(k, x1, u0, v)];
end
plot_x(x_tot, 'b');

if(compute_jacobians)
    tic
    Jacobian_u_G_temp = jacobian(G(lambda,u,x1,N,m,r,Q,c,Q_c,min_car_turning_radius,v), u);  % s x r
    Jacobian_u_J_temp = jacobian(J(u, N, x1, x_target, v), u); % 1 x r
    Rho_a_sym_temp = Rho(lambda, u, mu, a, N, r, s, x1, Jacobian_u_J_temp, Jacobian_u_G_temp, m, Q, c, Q_c,min_car_turning_radius,v); % (r + s) x 1
    
    % (r + s) x (r + s + 1)
    Jacobian_Rho_a_temp = jacobian(Rho(lambda, u, mu, a, N, r, s, x1, Jacobian_u_J_temp, Jacobian_u_G_temp, m, Q, c, Q_c, min_car_turning_radius,v), [lambda; u; mu]);
    Jacobian_u_G = Function('Jacobian_u_G',{lambda,u,mu},{Jacobian_u_G_temp},{'lambda','u','mu'},{'r'});
    Jacobian_u_J = Function('Jacobian_u_J',{lambda,u,mu},{Jacobian_u_J_temp},{'lambda','u','mu'},{'r'});
    Rho_a_sym = Function('Rho_a_sym',{lambda,u,mu},{Rho_a_sym_temp},{'lambda','u','mu'},{'r'});
    Jacobian_Rho_a = Function('Jacobian_Rho_a',{lambda,u,mu},{Jacobian_Rho_a_temp},{'lambda','u','mu'},{'r'});

    % Final step symbolic functions
    Rho_a_final_step = Function('Rho_a_final_step',{u,mu},{Rho(1, u, mu, a, N, r, s, x1, Jacobian_u_J_temp, Jacobian_u_G_temp, m, Q, c, Q_c, min_car_turning_radius,v)},{'u','mu'},{'r'});
    Jacobian_Rho_a_final_step = Function('Jacobian_Rho_a_final_step',{u,mu},{jacobian(Rho(1, u, mu, a, N, r, s, x1, Jacobian_u_J_temp, Jacobian_u_G_temp, m, Q, c, Q_c, min_car_turning_radius, v), [u; mu])},{'u','mu'},{'r'});
    time_to_get_symbolics = toc;
end

% find mu0 \in R^s_{>0}\
mu0 = getMu0(a, r, s, G(0,u0,x1,N,m,r,Q,c,Q_c,min_car_turning_radius,v));

if(verbose)
    temp_value = full(norm(Rho_a_sym(0, u0, mu0)));
    if(temp_value < 1e-6)
        disp('mu0 is good!');
    else
        disp('mu0 is not good!');
        disp(strcat('norm(Rho_a_sym(0, u0, mu0)) is :' , string(temp_value)));
    end
end

w_corrected = [0; u0; mu0]; % = [lambda; u; mu]
w_tot_corrected = w_corrected;
w_tot_predict = [];
direction_to_traverse_curve = 'unset';

tic
for i = 0:max_steps
    if(verbose)
        disp(strcat('///////////// i = ', string(i), '//////////////'));
    end

    [direction_to_traverse_curve, t] = getTangent_QR(full(evalf(Jacobian_Rho_a(w_corrected(1),  w_corrected(2:r+1), w_corrected(r+2:end)))), direction_to_traverse_curve);
    
    if(verbose)
        temp_value = full(norm(Jacobian_Rho_a(w_corrected(1),  w_corrected(2:r+1), w_corrected(r+2:end))*t));
        if(temp_value < 1e-6)
            disp('t is good!')
        else
            disp(strcat('t is not good! Norm is : ', string(norm(temp_value))));
        end
    end

    w_predict = w_corrected + ds*t;
    w_tot_predict = [w_tot_predict, w_predict];

    if((1 - full(w_predict(1))) < 0.01)
        tol = 1e-6;
        y_init = full(w_predict(2:end));
        y_final = Newton_corrector_for_final_step(y_init, Rho_a_final_step, Jacobian_Rho_a_final_step, r, max_iterations, tol);
        time_to_track_cure = toc;
        w_tot_corrected = [w_tot_corrected, [1; y_final]];

        break;
    else
        % Correct
        max_iterations = 100;
        tol = 1e-8;
        
        w_corrected = Newton_corrector(w_predict, Rho_a_sym, Jacobian_Rho_a, r, max_iterations, tol);
        
        if(verbose)
            temp_value = full(norm(Rho_a_sym(w_corrected(1), w_corrected(2:r+1), w_corrected(r+2:end))));
            if(temp_value < tol)
                disp('w_corrected is good!');
            else
                disp(strcat('w_corrected is not good! Norm is : ', string(temp_value)));
            end
            disp(strcat('lambda = ', string(full(w_corrected(1)))));
        end
        w_tot_corrected = [w_tot_corrected, w_corrected];
    end
end
% get the solution and plot
u_sol = w_tot_corrected(2:r+1,end);
mu_sol = w_tot_corrected(r+2:end,end);

x_tot = [];
for k=1:N+1
    x_tot = [x_tot, sol_diff_eqn(k, x1, u_sol, v)];
end

figure(1);
hold on;
plot_x(full(x_tot), color);
xlabel('$x_1$', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);
ylabel('$x_2$', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);
axis equal;
grid on;

% disp('Rho at end:');
% disp(Rho(1, u_sol, mu_sol, a, N, r, s, x1, Jacobian_u_J(1, u_sol, mu_sol), Jacobian_u_G(1, u_sol, mu_sol), m, Q, c, Q_c, min_car_turning_radius));
% 
% disp('G at end:');
% disp(G(1,u_sol,x1,N,m,r,Q,c, Q_c));

disp(strcat('Time to get symbolic expressions: ', string(time_to_get_symbolics)));

disp(strcat('Took: ', string(i), ' iterations to get lambda: ', string(w_tot_corrected(1,end))));

disp(strcat('Final cost homotopy: ', string(J(u_sol, N, x1, x_target, v))));

disp(strcat('Time to track zero curve: ', string(time_to_track_cure)));

plot_u(u_sol);