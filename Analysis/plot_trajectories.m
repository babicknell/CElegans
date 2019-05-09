function plot_trajectories(strain1, strain2)
% Plot tractories from two strain on one set of axes.
% Input: strings, eg. strain1 = 'CZ10175' and strain2 = 'QH6084'.

end_only = false; % change to 'true' for no axon

S1 = load(['../Results/', strain1, '.mat']);
S1 = S1.S;
S2 = load(['../Results/', strain2, '.mat']);
S2 = S2.S;

figure(1);
plot([0, 0], [0, 250],'b--','LineWidth',2);
hold on
plot([0, 0], [-50, 0],'b-','LineWidth',2)

for k = 1:length(S1)
    x = S1(k).Coords(:, 1);
    y = S1(k).Coords(:, 2);
    plot(x(end), y(end), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
    if ~end_only
        plot(x, y, 'k', 'LineWidth', 0.5)
    end
end

for k = 1:length(S2)
    x = S2(k).Coords(:, 1);
    y = S2(k).Coords(:, 2);
    plot(x(end), y(end), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r')
    if ~end_only
        plot(x, y, 'r', 'LineWidth', 0.5)
    end
end

axis equal
axis([-50, 50, -50, 200])
xlabel('displacement ($\mu{m}$)', 'Interpreter','Latex', 'FontSize',18)
ylabel({'original axon arclength ($\mu{m}$)'}, 'Interpreter', 'Latex', 'FontSize',18)
set(gca,'FontSize',20)