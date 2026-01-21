% Repl_Figs_LP.m
% -------------------------------------------------------------------------
% Local Projection version of the plotting code used for Figures


%% =========================
% Formatting (match the doc)
% =========================
fonttype      = 'Arial';
ftsizeaxis    = 11;
titlefontsize = 10;

% Map the "bound" (quantile tail prob used in est_irf_fev) to LP coverage (%)
% Example: bound=0.16 -> 68% bands
bandsCoverage = 100 * (1 - 2*bound);

%% =========================================================
% LP: set same variable ordering like replication
% =========================================================
choose_vars_LP = [12 7 2 3 4 5 9 10 8 11];

[data_lp, ~, vars, ~, ~, ~] = est_irf_lp( ...
    nirf, bound, ndraws, nconst, prior, begin_date, end_date, ...
    esty1, estq1, esty2, estq2, nlags, choose_vars_LP, samp, DATA_MASTER, labels_vars, ...
    'identification', 'CHOL', ...
    'shockVar', 1, ...
    'shockSize', 1, ...
    'bandsCoverage', bandsCoverage, ...
    'scale100', true);

plot_irf_grid_like(data_lp, vars, ...
    'Figure 6 (Local Projections)', fonttype, ftsizeaxis, titlefontsize);




% ============================
% Local function: plotting grid
% ============================
function plot_irf_grid_like(data, vars, figTitle, fonttype, ftsizeaxis, titlefontsize)

% data is (H+1) x N x 3: median/lb/ub
[H1, N, ~] = size(data);
H  = H1 - 1;

% --- Horizon axis: 0..H (recommended, consistent with XTick at 0,4,8,...)
hor = (0:H)';

R = round(N/2);

% Figure style like your original snippet
b = figure( ...
    'Color', [0.9412 0.9412 0.9412], ...
    'Position', [1 1 800-100 600-100], ...
    'Name', figTitle);
figure(b);

for i = 1:N
    subplot(R, 2, i);
    hold on;

    med = squeeze(data(:, i, 1));
    lb  = squeeze(data(:, i, 2));
    ub  = squeeze(data(:, i, 3));

    % --- Median (blue solid): do NOT set color; MATLAB default is blue
    plot(hor, med, 'LineWidth', 2);

    % --- Bands (red dashed)
    plot(hor, lb, 'r--', 'LineWidth', 2);
    plot(hor, ub, 'r--', 'LineWidth', 2);

    % --- Zero line
    plot(hor, zeros(size(hor)), ':k', 'LineWidth', 0.3);

    % Labels/titles
    ylabel('percent', 'FontSize', 12);
    xlabel('Horizon (quarters)', 'FontSize', 12);

    if ischar(vars)
        vname = strtrim(vars(i,:));
        title(vname, 'Interpreter','tex', 'FontSize', titlefontsize);
    else
        title(string(vars(i)), 'Interpreter','tex', 'FontSize', titlefontsize);
    end

    % X ticks like your original code, but robust if H < 20
    xt = 0:4:min(20,H);
    set(gca, 'XTick', xt);
    set(gca, 'XTickLabel', string(xt));

    % Font & cosmetics
    set(gca, 'FontName', fonttype);
    set(gca, 'FontSize', ftsizeaxis);
    set(gca, 'Layer', 'top');
    box off;
    axis tight;

    hold off;
end

end