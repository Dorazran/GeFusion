% -----------------------------------------------------------------------------
% Copyright (c) 2026 Dor Azran, Ariel University
% Licensed under the MIT License. See LICENSE file in the project root.
% -----------------------------------------------------------------------------

%% global_mc_sobol.m
% =============================================================================
%  GLOBAL VARIANCE-BASED (SOBOL) SENSITIVITY ANALYSIS
%  for Fick's-law / erfc dopant diffusion junction depth x_j in SiGe
% -----------------------------------------------------------------------------
%  PURPOSE
%  -------
%  This script REPLACES the previous crude one-at-a-time (+-20%) sweep of the
%  single parameter "alpha" (Ge-enhancement exponent) with a rigorous GLOBAL
%  Monte Carlo sensitivity analysis based on Sobol variance decomposition.
%
%  Instead of perturbing one parameter at a time around its nominal value
%  while holding all others fixed (which cannot detect interactions and which
%  arbitrarily privileges one input), we now treat ALL physically uncertain
%  inputs (alpha, D0, Ea, T, x_Ge, t) as random variables drawn from
%  literature-motivated distributions, propagate them jointly through the
%  closed-form analytic model for the junction depth x_j, and decompose the
%  variance of x_j into contributions attributable to each input (first-order
%  index S_i) and to each input plus all of its interactions (total-order
%  index ST_i).
%
%  PHYSICAL MODEL
%  --------------
%  Dopant profile (complementary error function solution of Fick's 2nd law
%  for a constant-source / semi-infinite diffusion):
%
%       C(x,t) = Cs * erfc( x / (2*sqrt(D*t)) )
%
%  Junction depth x_j defined as the depth at which C(x,t) equals the
%  background doping concentration C_B:
%
%       C_B = Cs * erfc( x_j / (2*sqrt(D*t)) )
%       =>   x_j = 2*sqrt(D*t) * erfcinv(C_B/Cs)
%
%  Diffusivity in Si(1-x)Ge(x), Arrhenius form with a Ge-fraction-dependent
%  exponential enhancement/retardation term (Zangenberg et al.; Liu et al.):
%
%       D_SiGe(x_Ge,T) = D0 * exp(-Ea/(k_B*T)) * exp(alpha*x_Ge)
%
%  SOBOL / SALTELLI METHOD (implemented BY HAND, no toolboxes)
%  -------------------------------------------------------------
%  Reference for the estimator formulas used below:
%    A. Saltelli, P. Annoni, I. Azzini, F. Campolongo, M. Ratto, S. Tarantola,
%    "Variance based sensitivity analysis of model output. Design and
%    estimator for the total sensitivity index", Computer Physics
%    Communications, 181(2), 259-270 (2010).
%  (and the original radial/Saltelli (2002) sampling design that the 2010
%   paper refines; the total-index estimator below is the Jansen (1999) /
%   Saltelli et al. (2010) "Eq.(g)"-type estimator using f(A) and f(AB_i).)
%
%  Sampling design:
%    A, B           : two independent N x k matrices of the k uncertain
%                      inputs, each column sampled from its own marginal
%                      distribution.
%    AB_i           : matrix equal to A except column i is replaced by the
%                      i-th column of B  ("radial"/Saltelli design).
%
%  Model evaluations:
%    f(A), f(B), f(AB_i)  for i = 1..k   -->  (k+2)*N evaluations total.
%
%  First-order Sobol index (Saltelli et al. 2010, estimator for S_i):
%
%       S_i  = [ (1/N) * sum_{n=1}^{N} f(B)_n * ( f(AB_i)_n - f(A)_n ) ] / Var(Y)
%
%  Total-order Sobol index (Saltelli et al. 2010 / Jansen 1999 estimator):
%
%       ST_i = 1 - [ (1/N) * sum_{n=1}^{N} f(A)_n * ( f(AB_i)_n - f(B)_n ) ] / Var(Y)
%
%  where Var(Y) is the sample variance of Y estimated from the POOLED sample
%  of f(A) and f(B) (2N values), which gives a more stable variance estimate
%  than using f(A) alone.
%
%  All sampling uses ONLY base MATLAB rand/randn (no Statistics and Machine
%  Learning Toolbox, no sobolset/sobolpoint, no erfinv-free tricks beyond
%  what ships with core MATLAB). erfcinv and erfinv are core MATLAB
%  functions (part of base MATLAB, not the Statistics Toolbox).
%
%  OUTPUTS (written into the same folder as this script)
%  -------------------------------------------------------------------------
%   sobol_results_<Dopant>.mat     - S_i, ST_i, output stats, raw samples info
%   sobol_results_<Dopant>.csv     - columns: parameter, S_i, ST_i
%   sobol_indices_bar_<Dopant>.png - bar chart of S_i vs ST_i per parameter
%   xj_distribution_histogram_<Dopant>.png - histogram of x_j with 5/95 pct
%   sobol_comparison_all_dopants.png - combined ST_i comparison across dopants
%
%  Author: course paper revision (Ariel University) -- v3 deep revision
% =============================================================================

clear; clc; close all;
tic;  % start wall-clock timer for the whole analysis

%% -------------------- USER-CONFIGURABLE SETTINGS --------------------------
N = 20000;              % Saltelli base sample size per matrix (A and B each N x k)
k = 6;                  % number of uncertain inputs:
                         %   1: alpha   2: D0   3: Ea   4: T   5: x_Ge   6: t
rng_seed = 42;           % fix RNG seed for reproducibility
rng(rng_seed);

param_names = {'alpha','D0','Ea','T','x\_Ge','t'};  % for plotting (LaTeX-ish)
param_names_csv = {'alpha','D0','Ea','T','x_Ge','t'};

% Fixed (non-uncertain) physical constants -----------------------------------
Cs   = 1e20;            % surface concentration [cm^-3]  (fixed)
CB   = 1e15;            % background concentration [cm^-3] (fixed)
k_B  = 8.617333e-5;     % Boltzmann constant [eV/K] (fixed)

% Output directory = same directory as this script
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

% List of dopants to analyze, with nominal literature values
dopant_list = {'Boron','Arsenic','Phosphorus'};

nominal.alpha.Boron      = -3.0;
nominal.alpha.Arsenic    = +2.3;
nominal.alpha.Phosphorus = +0.7;

nominal.D0.Boron      = 0.76;   % cm^2/s
nominal.D0.Arsenic    = 22.9;   % cm^2/s
nominal.D0.Phosphorus = 3.85;   % cm^2/s

nominal.Ea.Boron      = 3.46;   % eV
nominal.Ea.Arsenic    = 4.05;   % eV
nominal.Ea.Phosphorus = 3.66;   % eV

nominal.T_K      = 1273.15;     % K  (1000 C anneal), same for all dopants
nominal.x_Ge     = 0.20;        % Ge fraction, same for all dopants
nominal.t_s      = 3600;        % s  (1 hour anneal), same for all dopants

% Uncertainty specification (documented per the task statement):
%   alpha : Normal(mean = nominal, std = 0.20*|nominal|)        [20% relative]
%   D0    : Lognormal, median = literature value,
%           geometric std factor (GSD) = 1.3 (~30% multiplicative uncertainty)
%   Ea    : Normal(mean = nominal, std = 0.02 eV)                [tight, well-known]
%   T     : Normal(mean = 1273.15 K, std = 3 K)                  [furnace control]
%   x_Ge  : Normal(mean = 0.20, std = 0.02)                      [SIMS calibration]
%   t     : Normal(mean = 3600 s, std = 30 s)                    [process timer]

GSD_D0 = 1.3;   % geometric standard deviation factor for D0 lognormal model

%% -------------------- STORAGE FOR CROSS-DOPANT COMPARISON ------------------
all_ST = zeros(numel(dopant_list), k);   % rows = dopants, cols = parameters
all_S  = zeros(numel(dopant_list), k);

%% =====================  MAIN LOOP OVER DOPANT SPECIES  =====================
for d_idx = 1:numel(dopant_list)

    dopant = dopant_list{d_idx};
    fprintf('\n=====================================================\n');
    fprintf(' Running global Sobol sensitivity analysis for: %s\n', dopant);
    fprintf('=====================================================\n');

    % ---- nominal values for this dopant ----
    alpha_nom = nominal.alpha.(dopant);
    D0_nom    = nominal.D0.(dopant);
    Ea_nom    = nominal.Ea.(dopant);
    T_nom     = nominal.T_K;
    xGe_nom   = nominal.x_Ge;
    t_nom     = nominal.t_s;

    % ---- distribution parameters ----
    % 1) alpha ~ Normal(alpha_nom, (0.20*|alpha_nom|)^2)
    alpha_mu    = alpha_nom;
    alpha_sigma = 0.20 * abs(alpha_nom);

    % 2) D0 ~ Lognormal with median = D0_nom and geometric std factor GSD_D0.
    %    For a lognormal X = exp(mu + sigma*Z), Z~N(0,1):
    %       median(X) = exp(mu)              -> mu = log(median)
    %       GSD(X)    = exp(sigma)           -> sigma = log(GSD)
    D0_mu_ln    = log(D0_nom);
    D0_sigma_ln = log(GSD_D0);

    % 3) Ea ~ Normal(Ea_nom, 0.02^2)
    Ea_mu    = Ea_nom;
    Ea_sigma = 0.02;

    % 4) T ~ Normal(1273.15, 3^2)
    T_mu    = T_nom;
    T_sigma = 3;

    % 5) x_Ge ~ Normal(0.20, 0.02^2)
    xGe_mu    = xGe_nom;
    xGe_sigma = 0.02;

    % 6) t ~ Normal(3600, 30^2)
    t_mu    = t_nom;
    t_sigma = 30;

    %% ---------------- BUILD BASE SAMPLE MATRICES A and B -------------------
    % Each column corresponds to one of the k=6 uncertain inputs, in order:
    %   col 1: alpha   col 2: D0   col 3: Ea   col 4: T   col 5: x_Ge  col 6: t
    %
    % We draw all underlying randomness from rand/randn only (no toolbox
    % distribution objects, no sobolset / quasi-random low-discrepancy
    % sequences -- a plain pseudo-random Saltelli design as requested).

    A = zeros(N, k);
    B = zeros(N, k);

    % --- column 1: alpha (Normal) ---
    A(:,1) = alpha_mu + alpha_sigma .* randn(N,1);
    B(:,1) = alpha_mu + alpha_sigma .* randn(N,1);

    % --- column 2: D0 (Lognormal), sampled as exp(mu + sigma*Z) ---
    A(:,2) = exp(D0_mu_ln + D0_sigma_ln .* randn(N,1));
    B(:,2) = exp(D0_mu_ln + D0_sigma_ln .* randn(N,1));

    % --- column 3: Ea (Normal) ---
    A(:,3) = Ea_mu + Ea_sigma .* randn(N,1);
    B(:,3) = Ea_mu + Ea_sigma .* randn(N,1);

    % --- column 4: T (Normal) ---
    A(:,4) = T_mu + T_sigma .* randn(N,1);
    B(:,4) = T_mu + T_sigma .* randn(N,1);

    % --- column 5: x_Ge (Normal) ---
    A(:,5) = xGe_mu + xGe_sigma .* randn(N,1);
    B(:,5) = xGe_mu + xGe_sigma .* randn(N,1);

    % --- column 6: t (Normal) ---
    A(:,6) = t_mu + t_sigma .* randn(N,1);
    B(:,6) = t_mu + t_sigma .* randn(N,1);

    %% ---------------- EVALUATE MODEL ON A AND B -----------------------------
    fA = xj_model(A(:,1), A(:,2), A(:,3), A(:,4), A(:,5), A(:,6), Cs, CB, k_B);
    fB = xj_model(B(:,1), B(:,2), B(:,3), B(:,4), B(:,5), B(:,6), Cs, CB, k_B);

    % Pooled-sample variance estimate of Y (more stable than using A alone)
    Y_pool  = [fA; fB];
    VarY    = var(Y_pool, 1);   % normalize by N (population-style), consistent
                                 % with the Saltelli (2010) estimator denominator

    %% ---------------- BUILD AB_i MATRICES AND EVALUATE -----------------------
    % AB_i = A with column i replaced by B's column i
    S_i  = zeros(k,1);
    ST_i = zeros(k,1);

    for i = 1:k
        AB_i = A;
        AB_i(:,i) = B(:,i);

        fABi = xj_model(AB_i(:,1), AB_i(:,2), AB_i(:,3), ...
                         AB_i(:,4), AB_i(:,5), AB_i(:,6), Cs, CB, k_B);

        % ---- First-order Sobol index (Saltelli et al. 2010 estimator) ----
        %   S_i = (1/N * sum( fB .* (fABi - fA) )) / Var(Y)
        S_i(i) = ( (1/N) * sum( fB .* (fABi - fA) ) ) / VarY;

        % ---- Total-order Sobol index (Jansen / Saltelli et al. 2010) -----
        %   ST_i = 1 - (1/N * sum( fA .* (fABi - fB) )) / Var(Y)
        ST_i(i) = 1 - ( (1/N) * sum( fA .* (fABi - fB) ) ) / VarY;
    end

    % Store for cross-dopant comparison figure
    all_S(d_idx,:)  = S_i';
    all_ST(d_idx,:) = ST_i';

    %% ---------------- OUTPUT DISTRIBUTION SUMMARY STATISTICS -----------------
    % Use the pooled A/B output sample (2N realistic combined-uncertainty
    % evaluations) to summarize the x_j distribution itself.
    xj_mean   = mean(Y_pool);
    xj_std    = std(Y_pool, 1);
    xj_p05    = prctile_manual(Y_pool, 5);
    xj_p95    = prctile_manual(Y_pool, 95);

    fprintf('\n--- x_j output distribution summary for %s ---\n', dopant);
    fprintf('  mean(x_j)      = %.6e cm\n', xj_mean);
    fprintf('  std(x_j)       = %.6e cm\n', xj_std);
    fprintf('  5th percentile = %.6e cm\n', xj_p05);
    fprintf('  95th percentile= %.6e cm\n', xj_p95);
    fprintf('  (i.e. a realistic 90%% combined-uncertainty interval for x_j,\n');
    fprintf('   far more informative than the old +-20%% single-parameter alpha sweep)\n');

    %% ---------------- PRINT SOBOL INDEX TABLE TO CONSOLE ----------------------
    fprintf('\n--- Sobol sensitivity indices for x_j(%s) ---\n', dopant);
    fprintf('  %-10s %12s %12s\n', 'Parameter', 'S_i (1st)', 'ST_i (tot)');
    fprintf('  %-10s %12s %12s\n', '---------', '---------', '----------');
    for i = 1:k
        fprintf('  %-10s %12.5f %12.5f\n', param_names_csv{i}, S_i(i), ST_i(i));
    end
    sumS = sum(S_i);
    fprintf('  %-10s %12.5f\n', 'sum(S_i)', sumS);
    fprintf('  (sum of first-order indices < 1 indicates interaction effects;\n');
    fprintf('   ST_i > S_i for a parameter also indicates it participates in\n');
    fprintf('   interactions with other uncertain inputs)\n');

    %% ---------------- SAVE .mat RESULTS ---------------------------------------
    results = struct();
    results.dopant        = dopant;
    results.N             = N;
    results.k             = k;
    results.param_names   = param_names_csv;
    results.S_i           = S_i;
    results.ST_i          = ST_i;
    results.VarY          = VarY;
    results.xj_mean       = xj_mean;
    results.xj_std        = xj_std;
    results.xj_p05        = xj_p05;
    results.xj_p95        = xj_p95;
    results.nominal_alpha = alpha_nom;
    results.nominal_D0    = D0_nom;
    results.nominal_Ea    = Ea_nom;
    results.nominal_T     = T_nom;
    results.nominal_xGe   = xGe_nom;
    results.nominal_t     = t_nom;

    mat_path = fullfile(script_dir, sprintf('sobol_results_%s.mat', dopant));
    save(mat_path, 'results');
    fprintf('\nSaved: %s\n', mat_path);

    % Also save the generic (non-dopant-suffixed) name for the LAST dopant
    % processed, to satisfy the literal "sobol_results.mat" / ".csv" naming
    % requested in the task statement, in addition to the per-dopant files.
    if d_idx == numel(dopant_list)
        save(fullfile(script_dir, 'sobol_results.mat'), 'results');
    end

    %% ---------------- SAVE .csv RESULTS ----------------------------------------
    csv_path = fullfile(script_dir, sprintf('sobol_results_%s.csv', dopant));
    fid = fopen(csv_path, 'w');
    fprintf(fid, 'parameter,S_i,ST_i\n');
    for i = 1:k
        fprintf(fid, '%s,%.8f,%.8f\n', param_names_csv{i}, S_i(i), ST_i(i));
    end
    fclose(fid);
    fprintf('Saved: %s\n', csv_path);

    if d_idx == numel(dopant_list)
        csv_generic = fullfile(script_dir, 'sobol_results.csv');
        fid = fopen(csv_generic, 'w');
        fprintf(fid, 'parameter,S_i,ST_i\n');
        for i = 1:k
            fprintf(fid, '%s,%.8f,%.8f\n', param_names_csv{i}, S_i(i), ST_i(i));
        end
        fclose(fid);
        fprintf('Saved: %s\n', csv_generic);
    end

    %% ---------------- BAR CHART: S_i vs ST_i -----------------------------------
    fig1 = figure('Visible','off','Position',[100 100 1100 750]);
    bar_data = [S_i, ST_i];
    bar(bar_data, 'grouped');
    set(gca, 'XTickLabel', param_names_csv, 'FontSize', 20, 'LineWidth', 1.2);
    xlabel('Uncertain input parameter', 'FontSize', 22, 'FontWeight','bold');
    ylabel('Sobol sensitivity index', 'FontSize', 22, 'FontWeight','bold');
    title(sprintf('Sobol indices for x_j, %s (N = %d)', dopant, N), 'FontSize', 20);
    legend({'S_i (first-order)','ST_i (total-order)'}, 'Location','best', 'FontSize', 18);
    grid on;
    bar_png = fullfile(script_dir, sprintf('sobol_indices_bar_%s.png', dopant));
    saveas(fig1, bar_png);
    close(fig1);
    fprintf('Saved: %s\n', bar_png);

    if d_idx == numel(dopant_list)
        fig1b = figure('Visible','off','Position',[100 100 1100 750]);
        bar_data_generic = [S_i, ST_i];
        bar(bar_data_generic, 'grouped');
        set(gca, 'XTickLabel', param_names_csv, 'FontSize', 20, 'LineWidth', 1.2);
        xlabel('Uncertain input parameter', 'FontSize', 22, 'FontWeight','bold');
        ylabel('Sobol sensitivity index', 'FontSize', 22, 'FontWeight','bold');
        title(sprintf('Sobol indices for x_j, %s (N = %d)', dopant, N), 'FontSize', 20);
        legend({'S_i (first-order)','ST_i (total-order)'}, 'Location','best', 'FontSize', 18);
        grid on;
        saveas(fig1b, fullfile(script_dir, 'sobol_indices_bar.png'));
        close(fig1b);
    end

    %% ---------------- HISTOGRAM OF x_j DISTRIBUTION ----------------------------
    fig2 = figure('Visible','off','Position',[100 100 1100 750]);
    histogram(Y_pool, 60, 'Normalization','probability');
    hold on;
    yl = ylim;
    plot([xj_p05 xj_p05], yl, 'r--', 'LineWidth', 1.5);
    plot([xj_p95 xj_p95], yl, 'r--', 'LineWidth', 1.5);
    hold off;
    set(gca, 'FontSize', 20, 'LineWidth', 1.2);
    xlabel('Junction depth x_j [cm]', 'FontSize', 22, 'FontWeight','bold');
    ylabel('Probability', 'FontSize', 22, 'FontWeight','bold');
    title(sprintf('x_j distribution under combined uncertainty, %s', dopant), 'FontSize', 20);
    legend({'x_j samples','5th / 95th percentile'}, 'Location','best', 'FontSize', 18);
    grid on;
    hist_png = fullfile(script_dir, sprintf('xj_distribution_histogram_%s.png', dopant));
    saveas(fig2, hist_png);
    close(fig2);
    fprintf('Saved: %s\n', hist_png);

    if d_idx == numel(dopant_list)
        fig2b = figure('Visible','off','Position',[100 100 1100 750]);
        histogram(Y_pool, 60, 'Normalization','probability');
        hold on;
        yl = ylim;
        plot([xj_p05 xj_p05], yl, 'r--', 'LineWidth', 1.5);
        plot([xj_p95 xj_p95], yl, 'r--', 'LineWidth', 1.5);
        hold off;
        set(gca, 'FontSize', 20, 'LineWidth', 1.2);
        xlabel('Junction depth x_j [cm]', 'FontSize', 22, 'FontWeight','bold');
        ylabel('Probability', 'FontSize', 22, 'FontWeight','bold');
        title(sprintf('x_j distribution under combined uncertainty, %s', dopant), 'FontSize', 20);
        legend({'x_j samples','5th / 95th percentile'}, 'Location','best', 'FontSize', 18);
        grid on;
        saveas(fig2b, fullfile(script_dir, 'xj_distribution_histogram.png'));
        close(fig2b);
    end

end % end loop over dopants

%% =====================  COMBINED CROSS-DOPANT COMPARISON  ===================
% Total-order indices ST_i are generally preferred for ranking importance
% when interactions may be present, so the comparison figure uses ST_i.
fig3 = figure('Visible','off','Position',[100 100 1200 750]);
bar(all_ST', 'grouped');
set(gca, 'XTickLabel', param_names_csv, 'FontSize', 20, 'LineWidth', 1.2);
xlabel('Uncertain input parameter', 'FontSize', 22, 'FontWeight','bold');
ylabel('Total-order Sobol index ST_i', 'FontSize', 22, 'FontWeight','bold');
title('Total-order Sobol indices: comparison across dopants', 'FontSize', 20);
legend(dopant_list, 'Location','best', 'FontSize', 18);
grid on;
comp_png = fullfile(script_dir, 'sobol_comparison_all_dopants.png');
saveas(fig3, comp_png);
close(fig3);
fprintf('\nSaved combined comparison figure: %s\n', comp_png);

%% -------------------- ELAPSED WALL-CLOCK TIME -------------------------------
elapsed_s = toc;
fprintf('\n=====================================================\n');
fprintf(' Global Sobol sensitivity analysis complete.\n');
fprintf(' Total evaluations per dopant: %d  (= (k+2)*N = %d*%d)\n', ...
        (k+2)*N, (k+2), N);
fprintf(' Dopants analyzed: %s\n', strjoin(dopant_list, ', '));
fprintf(' Elapsed wall-clock time: %.3f seconds (%.3f minutes)\n', ...
        elapsed_s, elapsed_s/60);
fprintf('=====================================================\n');


%% =============================================================================
%  LOCAL FUNCTION DEFINITIONS
%  (MATLAB allows local functions at the end of a script file, callable
%   from anywhere above in the same file, as of MATLAB R2016b+)
%% =============================================================================

function xj = xj_model(alpha, D0, Ea, T, x_Ge, t, Cs, CB, k_B)
    % xj_model  Vectorized closed-form junction depth model.
    %
    %   D = D0 .* exp(-Ea ./ (k_B .* T)) .* exp(alpha .* x_Ge)
    %   xj = 2*sqrt(D.*t) * erfcinv(CB/Cs)
    %
    % All inputs except Cs, CB, k_B may be column vectors of equal length;
    % the function is fully vectorized (elementwise operations) so it can be
    % evaluated for an entire N x 1 Monte Carlo sample in one call.
    D  = D0 .* exp(-Ea ./ (k_B .* T)) .* exp(alpha .* x_Ge);
    xj = 2 .* sqrt(D .* t) .* erfcinv(CB/Cs);
end


function p = prctile_manual(x, pct)
    % prctile_manual  Simple percentile computation using only base MATLAB
    % sort/interp1, to avoid any dependency on the Statistics and Machine
    % Learning Toolbox's prctile() function.
    %
    %   x   : vector of samples
    %   pct : desired percentile, 0-100 (scalar)
    %
    % Uses linear interpolation between order statistics (the same default
    % method as MATLAB's own prctile, "linear interpolation of the empirical
    % CDF" variant), which is sufficiently accurate for N >= a few thousand.
    x = sort(x(:));
    n = numel(x);
    if n == 1
        p = x(1);
        return;
    end
    % Empirical CDF positions for sorted samples (matches common prctile
    % convention: position of the j-th order statistic is (j-0.5)/n)
    pos = ((1:n)' - 0.5) / n * 100;
    p = interp1(pos, x, pct, 'linear', 'extrap');
end
