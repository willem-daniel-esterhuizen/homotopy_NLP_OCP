function y_final = Newton_corrector_for_final_step(y_init, Rho_a_final_step, Jacobian_Rho_a_final_step, r, max_iterations, tol)

y_newton_loop = y_init;
for j = 1:max_iterations
    if(full(norm(Rho_a_final_step(y_newton_loop(1:r), y_newton_loop(r+1:end)))) < tol)
        break;
    else
        y_newton_loop = y_newton_loop - full(Jacobian_Rho_a_final_step(y_newton_loop(1:r), y_newton_loop(r+1:end)))\full(Rho_a_final_step(y_newton_loop(1:r), y_newton_loop(r+1:end)));
        continue;
    end
end
if ( (j == max_iterations) && (full(norm(Rho_a_final_step(y_newton_loop(1:r), y_newton_loop(r+1:end)))) > tol) )
    disp('final step did not converge');
else
    disp(strcat('final step converged in: ', string(j), ' steps' ));
end
y_final = y_newton_loop;
