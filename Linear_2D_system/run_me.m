% add path to Casadi
addpath('C:\Users\wiles\Downloads\casadi-3.6.7-windows64-matlab2018b');
clear all;

verbose = 1; % 1 for more info

import casadi.*

x1 = [0;0]; % initial state
x_target = [8;7];

figure_number = 1;
figure(figure_number);
xlabel('$x_1$', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);
ylabel('$x_2$', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);
grid on;
hold on;
plot(x1(1), x1(2), 'squarek');
plot(x_target(1), x_target(2), '^k');

% ellipsoidal constraints
Q1 = 0.4*eye(2);
Q = [Q1, Q1];
c = [[2;3], [7;5]]; % ellipsoids' centres

plot_ellipsoids(Q, c, figure_number);

% control's Euclidean norm is constrained:
control_bound = 0.5;
Q_c = eye(2);

n = 2; % state dimension
m = 2; % control dimension
N = 30; % horizon length
p = size(Q,2)/2; % number of state constraints
q = 1; % number of control constraints

r = N*m;
s = N*(p + q);

lambda = SX.sym('lambda',1);
u = SX.sym('u', r);
mu = SX.sym('mu', s);

% select a random a
u0 = rand(r,1);

% plot the initial guess
x_tot_init = [];
for k=1:N
    x_tot_init = [x_tot_init, sol_diff_eqn(k, x1, u0, m, r)];
end
plot_x(x_tot_init, 'bx-', figure_number);

b0 = get_b0(u0,x1,N,m,r,s,Q,c,Q_c,control_bound);

if(G(0,u0,x1,N,m,r,Q,c,Q_c,control_bound) < b0)
    if(verbose)
        disp('b0 is good!');
    end
end
c0 = ones(s,1);
a = [u0;b0;c0]; % should be (r + 2*s) x 1 

% Get all the symbolic expressions

tic

Jacobian_u_G_temp = jacobian(G(lambda,u,x1,N,m,r,Q,c,Q_c,control_bound), u);  % s x r
Jacobian_u_J_temp = jacobian(J(u, N, x1, x_target, m, r), u); % 1 x r
Rho_a_sym_temp = Rho(lambda, u, mu, a, N, r, s, x1, Jacobian_u_J_temp, Jacobian_u_G_temp, m, Q, c, Q_c,control_bound); % (r + s) x 1

% (r + s) x (r + s + 1)
Jacobian_Rho_a_temp = jacobian(Rho(lambda, u, mu, a, N, r, s, x1, Jacobian_u_J_temp, Jacobian_u_G_temp, m, Q, c, Q_c,control_bound), [lambda; u; mu]);

Jacobian_u_G = Function('Jacobian_u_G',{lambda,u,mu},{Jacobian_u_G_temp},{'lambda','u','mu'},{'r'});
Jacobian_u_J = Function('Jacobian_u_J',{lambda,u,mu},{Jacobian_u_J_temp},{'lambda','u','mu'},{'r'});
Rho_a_sym = Function('Rho_a_sym',{lambda,u,mu},{Rho_a_sym_temp},{'lambda','u','mu'},{'r'});
Jacobian_Rho_a = Function('Jacobian_Rho_a',{lambda,u,mu},{Jacobian_Rho_a_temp},{'lambda','u','mu'},{'r'});

% Symbolic functions for final Newton corrector
Rho_a_final_step = Function('Rho_a_final_step',{u,mu},{Rho(1, u, mu, a, N, r, s, x1, Jacobian_u_J_temp, Jacobian_u_G_temp, m, Q, c, Q_c, control_bound)},{'u','mu'},{'r'});
Jacobian_Rho_a_final_step = Function('Jacobian_Rho_a_final_step',{u,mu},{jacobian(Rho(1, u, mu, a, N, r, s, x1, Jacobian_u_J_temp, Jacobian_u_G_temp, m, Q, c, Q_c, control_bound), [u; mu])},{'u','mu'},{'r'});

time_to_compute_syms = toc;

% find mu0 \in R^s_{>0}\
mu0 = getMu0(a, r, s, G(0,u0,x1,N,m,r,Q,c,Q_c,control_bound));

if(verbose)
    temp_value = full(norm(Rho_a_sym(0, u0, mu0)));
    if(temp_value < 1e-6)
        disp('mu0 is good!');
    else
        disp('mu0 is not good!');
        % disp(strcat('norm(Rho_a_sym(0, u0, mu0)) is :' , string(temp_value)));
    end
end

max_steps = 100;
ds = 0.5; % step size

w_corrected = [0; u0; mu0]; % = [lambda; u; mu]
w_tot_corrected = w_corrected;
w_tot_predict = [];

tic
for i = 0:max_steps
    if(verbose)
        disp(strcat('///////////// i = ', string(i), '//////////////'));
    end

    t = getTangent_QR(full(evalf(Jacobian_Rho_a(w_corrected(1),  w_corrected(2:r+1), w_corrected(r+2:end)))));
    
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

    if( (1 - full(w_predict(1)) < 0.01) || (i == max_steps) )
        % final corrector
        disp(strcat('i at end: ', string(i)));
        tol = 1e-6;
        y_init = full(w_predict(2:end));
        [convergence_flag, y_final] = Newton_corrector_for_final_step(y_init, Rho_a_final_step, Jacobian_Rho_a_final_step, r, max_iterations, tol);
    
        w_tot_corrected = [w_tot_corrected, [1; y_final]];
        time_to_track = toc;
        break;
    else
        % Correct
        max_iterations = 100;
        tol = 1e-3;
        
        w_corrected = Newton_corrector(w_predict, Rho_a_sym, Jacobian_Rho_a, r, max_iterations, tol);
        
        if(verbose)
            temp_value = full(norm(Rho_a_sym(w_corrected(1), w_corrected(2:r+1), w_corrected(r+2:end))));
            if(temp_value < 1e-3)
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
    x_tot = [x_tot, sol_diff_eqn(k, x1, u_sol, m, r)];
end
plot_x(full(x_tot), 'ro-', figure_number);

% disp('Rho at end:');
% disp(Rho(1, u_sol, mu_sol, a, N, r, s, x1, Jacobian_u_J(1, u_sol, mu_sol), Jacobian_u_G(1, u_sol, mu_sol), m, Q, c, Q_c));
% 
% disp('G at end:');
% disp(G(1,u_sol,x1,N,m,r,Q,c, Q_c));

% disp('lambda:')
% w_tot_corrected(1,:)

figure(2);
hold on;
plot(w_tot_corrected(1,:), 'ko-');
xlabel('iterations', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);
ylabel('$\lambda$', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);
grid on;

u_parsed = parse_u(u_sol, m, r);

disp(strcat('Time to set up symbolic expressions: ', string(time_to_compute_syms)));
disp(strcat('Time to track zero curve: ', string(time_to_track)));


