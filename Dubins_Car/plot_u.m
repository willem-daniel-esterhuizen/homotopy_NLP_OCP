function plot_u(u_sol)

figure(2);
hold on;
grid on;

for k = 1:length(u_sol)
    plot([(k-1)*0.1, (k)*0.1], [u_sol(k), u_sol(k)], 'k');
    if k == length(u_sol)
        break;
    else
        plot([k*0.1, k*0.1], [u_sol(k), u_sol(k+1)], 'k');
    end
end

xlabel('$time (seconds)$', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);
ylabel('$u$', 'interpreter','latex','FontName','Times New Roman', 'FontSize', 15);