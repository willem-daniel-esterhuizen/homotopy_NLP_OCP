function w_corrected = Newton_corrector(w_predict, Rho_a_sym, Jacobian_Rho_a, r, max_iterations, tol)

w_newton_loop = w_predict;
for j = 1:max_iterations
    if(full(norm(Rho_a_sym(w_newton_loop(1), w_newton_loop(2:r+1), w_newton_loop(r+2:end)))) < tol)
        break;
    else
        w_newton_loop = w_newton_loop - pinv(full(Jacobian_Rho_a(w_newton_loop(1), w_newton_loop(2:r+1), w_newton_loop(r+2:end))) ) * full(Rho_a_sym(w_newton_loop(1), w_newton_loop(2:r+1), w_newton_loop(r+2:end)));
        continue;
    end
end
w_corrected = w_newton_loop;

