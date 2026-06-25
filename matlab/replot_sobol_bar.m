% -----------------------------------------------------------------------------
% Copyright (c) 2026 Dor Azran, Ariel University
% Licensed under the MIT License. See LICENSE file in the project root.
% -----------------------------------------------------------------------------

%% replot_sobol_bar.m
% Regenerates sobol_indices_bar.png (Fig. 5 in the AESMT manuscript) from the
% already-computed Phosphorus Sobol results, with much larger/bolder axis
% text so the chart stays legible after being shrunk into the manuscript's
% narrow two-column body width.
script_dir = fileparts(mfilename('fullpath'));

S = load(fullfile(script_dir, 'sobol_results_Phosphorus.mat'));
results = S.results;

param_names_csv = results.param_names;
S_i  = results.S_i;
ST_i = results.ST_i;
dopant = 'Phosphorus';
N = results.N;

fig1 = figure('Visible','off','Position',[100 100 1300 950]);
bar_data = [S_i, ST_i];
bar(bar_data, 'grouped');
set(gca, 'XTickLabel', param_names_csv, 'FontSize', 28, 'FontWeight','bold', 'LineWidth', 1.6);
xlabel('Uncertain input parameter', 'FontSize', 30, 'FontWeight','bold');
ylabel('Sobol sensitivity index', 'FontSize', 30, 'FontWeight','bold');
title(sprintf('Sobol indices for x_j, %s (N = %d)', dopant, N), 'FontSize', 26, 'FontWeight','bold');
lg = legend({'S_i (first-order)','ST_i (total-order)'}, 'Location','best', 'FontSize', 24);
set(lg, 'Box', 'on');
grid on;
set(gca, 'GridAlpha', 0.4);

bar_png = fullfile(script_dir, 'sobol_indices_bar.png');
saveas(fig1, bar_png);
close(fig1);
fprintf('Saved: %s\n', bar_png);
