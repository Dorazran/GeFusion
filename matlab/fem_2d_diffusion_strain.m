% -----------------------------------------------------------------------------
% Copyright (c) 2026 Dor Azran, Ariel University
% Licensed under the MIT License. See LICENSE file in the project root.
% -----------------------------------------------------------------------------

%% fem_2d_diffusion_strain.m
% =========================================================================
%  2D STRAIN-COUPLED DOPANT DIFFUSION ON A FinFET-LIKE CROSS-SECTION
%  (Ariel University course paper -- SiGe nano-device diffusion, revision v3)
% =========================================================================
%
%  PHYSICAL MOTIVATION
%  --------------------
%  The rest of this paper's simulation suite treats dopant diffusion in a
%  Ge-graded SiGe layer as a 1D problem: the Ge mole fraction x_Ge varies
%  only with depth z, and Fick's second law is solved along that single
%  depth coordinate. That is a reasonable first approximation for a planar
%  blanket SiGe film, but real nano-scale devices (e.g. FinFETs) have a
%  Ge-graded region that is bounded laterally by a fin sidewall, and the
%  Ge profile (and the resulting lattice-mismatch strain field) is known
%  from process integration studies to be non-uniform near fin corners
%  and sidewalls -- the Ge condensation / grading process does not produce
%  a perfectly depth-only profile once the fin width becomes comparable to
%  the grading length scale.
%
%  This script extends the 1D treatment to a 2D (x = across-fin, y = into
%  substrate) finite-difference solve on a rectangular cross-section that
%  is meant to stand in for a simplified FinFET fin + underlying graded
%  SiGe region. The purpose is illustrative: to show, within a course-paper
%  scope, that a laterally-resolved geometry with a coupled strain field
%  can deviate meaningfully from the simplified constant-Ge / depth-only
%  1D approximation used elsewhere in the paper -- NOT to provide a
%  process-calibrated TCAD-grade prediction.
%
%  NUMERICAL METHOD (IMPORTANT -- READ BEFORE USE)
%  --------------------------------------------------
%  Despite the filename "fem_2d_...", this is implemented with a
%  finite-DIFFERENCE method on a regular structured grid, NOT a true
%  finite-element method, and NOT MATLAB's PDE Toolbox (which may not be
%  licensed on the grading/target machine). The file name is kept for
%  continuity with the rest of the project's naming/file-tracking scheme;
%  the method actually used is a 2D Alternating-Direction-Implicit (ADI)
%  finite-difference scheme, built only from core MATLAB (sparse matrices,
%  spdiags, and the backslash operator), which are part of base MATLAB and
%  do not require any toolbox.
%
%  Governing equation (Fick's 2nd law with spatially varying diffusivity):
%       dC/dt = d/dx( D_eff(x,y) dC/dx ) + d/dy( D_eff(x,y) dC/dy )
%
%  D_eff(x,y) depends on position through (a) the local Ge mole fraction
%  x_Ge(x,y), which modulates the diffusivity via the usual SiGe
%  diffusivity-enhancement exponential exp(alpha*x_Ge), and (b) a
%  phenomenological strain-coupling correction (1 + beta*eps(x,y)) derived
%  from the Vegard's-law lattice-mismatch strain (see below).
%
%  Each full time step is split into two implicit half-steps (ADI):
%    1) Implicit sweep in x: for each of the NY grid rows, solve a
%       tridiagonal linear system (Thomas-algorithm-equivalent, done here
%       via MATLAB's sparse \ on a tridiagonal sparse matrix) advancing
%       C by dt/2 using x-direction diffusion only, y-direction explicit.
%    2) Implicit sweep in y: for each of the NX grid columns, solve a
%       tridiagonal linear system advancing C by dt/2 using y-direction
%       diffusion only, x-direction explicit (using the half-step result).
%  This ADI/backward-Euler splitting is unconditionally stable, which is
%  why we are free to choose NT based on accuracy considerations alone
%  and do NOT need to satisfy the explicit-scheme CFL/stability limit
%  dt <= dx^2/(2D) that would otherwise bound the largest usable time step.
%
%  SIMPLIFICATIONS / ILLUSTRATIVE ELEMENTS (explicitly flagged)
%  ----------------------------------------------------------------
%  (1) The lateral "corner rounding" term in x_Ge(x,y) (the
%      exp(-((x-W/2)^2)/(2*sigma_corner^2)) factor) is a simplified,
%      illustrative way to capture the kind of corner/sidewall
%      non-uniformity typical of real FinFET process integration. It is
%      NOT extracted from a measured Ge profile or a calibrated process
%      simulation; it is included purely to make the 2D problem
%      non-trivial in the lateral direction and to demonstrate the
%      qualitative effect of lateral Ge grading on the diffusion field.
%  (2) The strain-diffusivity coupling coefficient beta is an
%      illustrative, order-of-magnitude phenomenological constant, NOT a
%      fitted or literature-derived value. Rigorous ab-initio or
%      experimentally calibrated strain-diffusivity coupling constants
%      for B/As/P in strained SiGe are beyond the scope of this course
%      paper; beta is included only to demonstrate, structurally, how a
%      strain field arising from Vegard's-law lattice mismatch could be
%      folded into the diffusivity, and to provide a qualitative sense of
%      the magnitude of the effect under a simple linear-correction model.
%
% =========================================================================

clear; clc; close all;
script_start_time = tic;

%% -------------------- USER-CONFIGURABLE PARAMETERS --------------------

% --- Dopant selection ---------------------------------------------------
dopant = 'Boron';   % 'Boron' (only Boron literature values supplied here;
                     % extend dopant_table below to add As/P consistent
                     % with the rest of the project if needed)

% --- Geometry (FinFET-like cross-section) --------------------------------
W   = 0.2e-4;   % domain width in x (across the fin), cm  (0.2 um)
Dp  = 0.3e-4;   % domain depth in y (into the substrate), cm (0.3 um)
NX  = 80;       % number of grid nodes in x
NY  = 120;      % number of grid nodes in y

% --- Ge grading profile parameters (same as existing 1D GE_GRADE_PARAMS) -
x0_Ge       = 0.05;     % Ge fraction at the surface / top of grade
x1_Ge       = 0.25;     % Ge fraction at the bottom of the graded region
W_base      = 0.10e-4;  % characteristic grading depth scale, cm (0.10 um)
sigma_corner= 0.03e-4;  % lateral corner-rounding length scale, cm (0.03 um)
                         % (illustrative -- see header comment, item (1))
corner_amp  = 0.3;       % amplitude of the corner-rounding reduction term

% --- Anneal conditions ----------------------------------------------------
T_anneal   = 1273.15;   % nominal anneal temperature, K (1000 C)
t_anneal   = 3600;      % total anneal time, s (1 hour)
NT         = 2000;      % number of ADI time steps (see header: ADI removes
                         % the explicit CFL restriction, so NT is chosen
                         % for temporal accuracy, not stability)

% --- Boundary / initial conditions ----------------------------------------
Cs = 1e20;   % cm^-3, constant surface concentration at y = 0
CB = 1e15;   % cm^-3, background (initial) concentration elsewhere

% --- Strain coupling -------------------------------------------------------
beta = 5;    % illustrative strain-diffusivity coupling coefficient
             % (order-of-magnitude only -- see header comment, item (2))

% --- Lattice constants for Vegard's law (cm) -------------------------------
a_Si = 5.431e-8;   % cm
a_Ge = 5.658e-8;   % cm

% --- Physical constants -----------------------------------------------------
k_B = 8.617333262e-5;  % eV/K (Boltzmann constant)

%% -------------------- DOPANT DIFFUSIVITY LITERATURE TABLE -----------------
% D0 (cm^2/s), Ea (eV), alpha (dimensionless Ge-enhancement exponent)
% Values consistent with the rest of the project's 1D solver inputs.
dopant_table = struct( ...
    'Boron', struct('D0', 0.76, 'Ea', 3.46, 'alpha', -3.0) );

if ~isfield(dopant_table, dopant)
    error('fem_2d_diffusion_strain:UnknownDopant', ...
        'Dopant "%s" not found in dopant_table. Add its D0/Ea/alpha values.', dopant);
end
D0    = dopant_table.(dopant).D0;
Ea    = dopant_table.(dopant).Ea;
alpha = dopant_table.(dopant).alpha;

fprintf('=== 2D Strain-Coupled Diffusion Solver (FD/ADI) ===\n');
fprintf('Dopant: %s   D0=%.4g cm^2/s   Ea=%.4g eV   alpha=%.4g\n', dopant, D0, Ea, alpha);
fprintf('Anneal: T=%.2f K (%.1f C), t=%.1f s, NT=%d steps\n', ...
    T_anneal, T_anneal-273.15, t_anneal, NT);
fprintf('Grid: NX=%d x NY=%d nodes,  W=%.3f um, Dp=%.3f um\n', ...
    NX, NY, W*1e4, Dp*1e4);

%% -------------------- GRID CONSTRUCTION ------------------------------------
x = linspace(0, W,  NX);   % cm
y = linspace(0, Dp, NY);   % cm
dx = x(2) - x(1);
dy = y(2) - y(1);
dt = t_anneal / NT;

[Xg, Yg] = meshgrid(x, y);   % Xg, Yg are NY x NX (rows=y, cols=x), MATLAB convention

%% -------------------- Ge FRACTION FIELD x_Ge(x,y) --------------------------
% x_Ge(x,y) = x0 + (x1-x0) * clamp(y/W_base, 0, 1) * (1 - corner_amp * ...
%                exp(-((x-W/2)^2)/(2*sigma_corner^2)))
%
% The depth-dependence (clamp(y/W_base,0,1)) reproduces the same Ge grading
% law used in the existing 1D solver (ramping from x0 at the surface to x1
% over a depth scale W_base, then saturating at x1 beyond that depth).
%
% The lateral factor (1 - corner_amp*exp(-((x-W/2)^2)/(2*sigma_corner^2)))
% is the SIMPLIFIED, ILLUSTRATIVE corner-rounding term described in the
% header: it locally suppresses the Ge grading near the fin center-line
% (x = W/2 here is used purely as a stand-in coordinate for "near the
% sidewall-affected region"; the exact placement is arbitrary/illustrative)
% to represent the kind of non-uniform corner effect seen qualitatively in
% real FinFET Ge condensation/grading processes. It is NOT a measured or
% process-calibrated profile.

depth_ramp = min(max(Yg / W_base, 0), 1);                 % NY x NX, in [0,1]
corner_term = 1 - corner_amp * exp(-((Xg - W/2).^2) / (2*sigma_corner^2));
x_Ge = x0_Ge + (x1_Ge - x0_Ge) .* depth_ramp .* corner_term;   % NY x NX

% Average Ge fraction over the whole domain, used later for the 1D
% constant-Ge reference comparison.
x_Ge_avg = mean(x_Ge(:));
fprintf('Average Ge fraction over domain: %.4f\n', x_Ge_avg);

%% -------------------- STRAIN FIELD VIA VEGARD'S LAW -------------------------
% Local relaxed SiGe lattice constant via Vegard's law (linear interpolation
% between a_Si and a_Ge), and the resulting lattice-mismatch strain relative
% to unstrained Si:
%   a_SiGe(x,y) = a_Si + x_Ge(x,y) * (a_Ge - a_Si)
%   eps(x,y)    = (a_SiGe(x,y) - a_Si) / a_Si  =  x_Ge(x,y) * (a_Ge-a_Si)/a_Si
a_SiGe = a_Si + x_Ge .* (a_Ge - a_Si);
eps_strain = (a_SiGe - a_Si) ./ a_Si;     % NY x NX, dimensionless strain field

%% -------------------- DIFFUSIVITY FIELD D_eff(x,y) ---------------------------
% Base (Ge-fraction-dependent) Arrhenius diffusivity:
%   D(x,y,T) = D0 * exp(-Ea/(k_B*T)) * exp(alpha * x_Ge(x,y))
% Strain-corrected effective diffusivity (phenomenological linear coupling,
% see header item (2)):
%   D_eff(x,y) = D(x,y,T) * (1 + beta*eps(x,y))
D_arrhenius_prefactor = D0 * exp(-Ea / (k_B * T_anneal));   % scalar, cm^2/s
D_field   = D_arrhenius_prefactor * exp(alpha * x_Ge);       % NY x NX, cm^2/s
strain_correction = (1 + beta * eps_strain);

% Guard against an unphysical sign flip of the strain correction factor
% (possible only for extreme/unphysical beta*eps combinations); clip at a
% small positive floor so the diffusivity field stays well-posed.
strain_correction = max(strain_correction, 0.05);

D_eff = D_field .* strain_correction;   % NY x NX, cm^2/s -- the field actually used

fprintf('D_eff range over domain: [%.4g, %.4g] cm^2/s\n', min(D_eff(:)), max(D_eff(:)));

%% -------------------- INITIAL CONDITION --------------------------------------
% C(x,y,0) = CB everywhere except the y=0 surface row, which is held at Cs
% (Dirichlet boundary, enforced every step below as well).
C = CB * ones(NY, NX);
C(1, :) = Cs;   % y=0 is row index 1 in this (row=y, col=x) convention

%% -------------------- PRECOMPUTE ADI HALF-STEP COEFFICIENT ARRAYS ------------
% We use harmonic-mean face diffusivities for the variable-coefficient
% finite-difference Laplacian, which is the standard robust choice for
% spatially varying diffusivity on a structured grid.
%
% For a 1D term d/dx( D dC/dx ) discretized at node i (interior), with
% face diffusivities D_{i+1/2} and D_{i-1/2}:
%   d/dx(D dC/dx) |_i  ~=  [ D_{i+1/2}(C_{i+1}-C_i) - D_{i-1/2}(C_i-C_{i-1}) ] / dx^2
%
% Harmonic mean of D at a face avoids overestimating flux across regions of
% rapidly varying diffusivity (standard FD practice for variable-coefficient
% diffusion problems).
harm_mean = @(Da, Db) 2*Da.*Db ./ (Da + Db + eps);

%% -------------------- TIME-MARCHING LOOP (ADI / BACKWARD EULER) --------------
% Note on stability: because this is a fully-implicit ADI scheme, it is
% unconditionally stable in the usual von Neumann sense for the linear
% diffusion operator, regardless of how large dt is relative to dx^2/D.
% This is precisely why we can choose NT based on desired temporal
% resolution/accuracy rather than the strict CFL-type limit dt <= dx^2/(2D)
% that would be required for a fully explicit FTCS scheme.

progress_marks = round((1:10) * NT / 10);

for n = 1:NT

    % =====================================================================
    % HALF-STEP 1: implicit in x, explicit in y  (advance by dt/2)
    % =====================================================================
    C_half = zeros(NY, NX);

    for j = 1:NY   % loop over rows (fixed y), solve a tridiagonal system in x

        % Build the explicit y-direction contribution for this row using
        % values at rows j-1, j, j+1 (Neumann / zero-flux handled via
        % ghost-node mirroring at j=1 and j=NY).
        if j == 1
            % y=0 row is a Dirichlet BC row (Cs); no PDE update needed here,
            % but we still need a valid C_half row to feed into the y-sweep.
            C_half(j, :) = Cs;
            continue;
        end

        jp = min(j+1, NY);
        jm = max(j-1, 1);

        if j == NY
            % zero-flux (Neumann) at y = Dp: mirror ghost node, i.e. use
            % the same diffusivity/concentration as the interior neighbor
            % so that the flux term in y vanishes for this explicit part.
            Dy_p = D_eff(j, :);          % effectively zero contribution below
            Dy_m = harm_mean(D_eff(j,:), D_eff(jm,:));
            Cy_term = Dy_m .* (C(jm, :) - C(j, :)) / dy^2;
        else
            Dy_p = harm_mean(D_eff(j,:), D_eff(jp,:));
            Dy_m = harm_mean(D_eff(j,:), D_eff(jm,:));
            Cy_term = ( Dy_p .* (C(jp,:) - C(j,:)) - Dy_m .* (C(j,:) - C(jm,:)) ) / dy^2;
        end

        % Right-hand side for this row's implicit x-tridiagonal solve:
        % C_row^{n+1/2} - (dt/2)*d/dx(D dC/dx)|^{n+1/2} = C_row^n + (dt/2)*Cy_term^n
        rhs = C(j, :)' + (dt/2) * Cy_term';   % NX x 1

        % Build tridiagonal matrix for implicit x-direction operator on this row.
        Drow = D_eff(j, :);  % 1 x NX diffusivity along this row

        main_diag = ones(NX,1);
        lower     = zeros(NX,1);
        upper     = zeros(NX,1);

        for i = 2:NX-1
            Dxp = 2*Drow(i)*Drow(i+1) / (Drow(i)+Drow(i+1)+eps);
            Dxm = 2*Drow(i)*Drow(i-1) / (Drow(i)+Drow(i-1)+eps);
            cE  = (dt/2) * Dxp / dx^2;
            cW  = (dt/2) * Dxm / dx^2;
            lower(i)     = -cW;
            main_diag(i) = 1 + cW + cE;
            upper(i)     = -cE;
        end

        % Neumann (zero-flux) BC at x=0 (i=1): mirror node => only flux to i=2.
        Dxp1 = 2*Drow(1)*Drow(2) / (Drow(1)+Drow(2)+eps);
        cE1  = (dt/2) * Dxp1 / dx^2;
        main_diag(1) = 1 + cE1;
        upper(1)     = -cE1;

        % Neumann (zero-flux) BC at x=W (i=NX): mirror node => only flux to i=NX-1.
        DxmN = 2*Drow(NX)*Drow(NX-1) / (Drow(NX)+Drow(NX-1)+eps);
        cWN  = (dt/2) * DxmN / dx^2;
        main_diag(NX) = 1 + cWN;
        lower(NX)     = -cWN;

        A_x = spdiags([ [lower(2:NX); 0], main_diag, [0; upper(1:NX-1)] ], [-1, 0, 1], NX, NX);
        % NOTE: spdiags expects each column already aligned to the diagonal
        % offset convention; lower(2:NX) holds subdiagonal entries (rows 2..NX),
        % upper(1:NX-1) holds superdiagonal entries (rows 1..NX-1). Padding
        % with a trailing/leading 0 keeps each column length = NX as required.

        C_half(j, :) = (A_x \ rhs)';
    end

    % =====================================================================
    % HALF-STEP 2: implicit in y, explicit in x  (advance remaining dt/2)
    % =====================================================================
    C_new = zeros(NY, NX);
    C_new(1, :) = Cs;   % enforce Dirichlet surface BC at y=0

    for i = 1:NX   % loop over columns (fixed x), solve a tridiagonal system in y

        ip = min(i+1, NX);
        im = max(i-1, 1);

        if i == 1
            Dx_p = harm_mean(D_eff(:,i), D_eff(:,ip));
            Cx_term = Dx_p .* (C_half(:,ip) - C_half(:,i)) / dx^2;   % zero-flux at i=1
        elseif i == NX
            Dx_m = harm_mean(D_eff(:,i), D_eff(:,im));
            Cx_term = -Dx_m .* (C_half(:,i) - C_half(:,im)) / dx^2;  % zero-flux at i=NX
        else
            Dx_p = harm_mean(D_eff(:,i), D_eff(:,ip));
            Dx_m = harm_mean(D_eff(:,i), D_eff(:,im));
            Cx_term = ( Dx_p .* (C_half(:,ip) - C_half(:,i)) - Dx_m .* (C_half(:,i) - C_half(:,im)) ) / dx^2;
        end

        rhs = C_half(:, i) + (dt/2) * Cx_term;   % NY x 1

        Dcol = D_eff(:, i);  % NY x 1 diffusivity along this column

        main_diag = ones(NY,1);
        lower     = zeros(NY,1);
        upper     = zeros(NY,1);

        for j = 2:NY-1
            Dyp = 2*Dcol(j)*Dcol(j+1) / (Dcol(j)+Dcol(j+1)+eps);
            Dym = 2*Dcol(j)*Dcol(j-1) / (Dcol(j)+Dcol(j-1)+eps);
            cN  = (dt/2) * Dyp / dy^2;   % toward larger y (deeper)
            cS  = (dt/2) * Dym / dy^2;   % toward smaller y (shallower)
            lower(j)     = -cS;
            main_diag(j) = 1 + cS + cN;
            upper(j)     = -cN;
        end

        % Dirichlet BC at y=0 (j=1): C fixed at Cs.
        main_diag(1) = 1;
        upper(1)     = 0;
        rhs(1)       = Cs;

        % Neumann (zero-flux) BC at y=Dp (j=NY): mirror node => only flux to j=NY-1.
        DymN = 2*Dcol(NY)*Dcol(NY-1) / (Dcol(NY)+Dcol(NY-1)+eps);
        cSN  = (dt/2) * DymN / dy^2;
        main_diag(NY) = 1 + cSN;
        lower(NY)     = -cSN;

        A_y = spdiags([ [lower(2:NY); 0], main_diag, [0; upper(1:NY-1)] ], [-1, 0, 1], NY, NY);

        C_new(:, i) = A_y \ rhs;
    end

    C = C_new;
    C(1, :) = Cs;   % re-assert Dirichlet BC exactly (guards against round-off)

    if any(n == progress_marks)
        pct = 100 * n / NT;
        fprintf('  progress: %3.0f%%  (step %d / %d, t = %.1f s)\n', pct, n, NT, n*dt);
    end
end

fprintf('Time-marching complete.\n');

%% -------------------- 1D CONSTANT-Ge REFERENCE PROFILE (erfc) ----------------
% Build a simple constant-Ge 1D reference using the AVERAGE Ge fraction
% over the whole 2D domain, to show how much the 2D / laterally-graded
% treatment deviates from the simplified constant-Ge approximation that
% would result from collapsing the lateral variation away.
D_ref_prefactor = D0 * exp(-Ea / (k_B * T_anneal));
D_ref = D_ref_prefactor * exp(alpha * x_Ge_avg);   % scalar effective diffusivity, cm^2/s

% Classic constant-source-concentration erfc diffusion profile:
%   C(y,t) = CB + (Cs - CB) * erfc( y / (2*sqrt(D_ref*t)) )
C_1d_ref = CB + (Cs - CB) * erfc( y / (2*sqrt(D_ref * t_anneal)) );   % 1 x NY

%% -------------------- CENTERLINE EXTRACTION FROM 2D SOLUTION -----------------
i_center = round(NX/2);
C_centerline = C(:, i_center)';   % 1 x NY, along x = W/2

%% -------------------- RELATIVE DIFFERENCE METRIC ------------------------------
% Avoid division-by-zero / domination by background noise: compute relative
% percentage difference only where the reference concentration is
% meaningfully above background (here, above 10*CB).
mask = C_1d_ref > 10*CB;
rel_diff_pct = 100 * abs(C_centerline(mask) - C_1d_ref(mask)) ./ C_1d_ref(mask);
max_rel_diff_pct = max(rel_diff_pct);

fprintf('\n=== 2D (laterally-graded, strain-coupled) vs 1D constant-Ge reference ===\n');
fprintf('Average Ge fraction used for 1D reference: %.4f\n', x_Ge_avg);
fprintf('Max relative percentage difference (centerline vs 1D erfc reference): %.2f %%\n', max_rel_diff_pct);

%% -------------------- SAVE RESULTS --------------------------------------------
out_dir = fileparts(mfilename('fullpath'));
if isempty(out_dir)
    out_dir = pwd;
end

results_file = fullfile(out_dir, 'fem2d_results.mat');
save(results_file, 'x', 'y', 'C', 'D_eff', 'x_Ge', 'eps_strain', ...
     'C_centerline', 'C_1d_ref', 'x_Ge_avg', 'D_ref', 'max_rel_diff_pct', ...
     'dopant', 'T_anneal', 't_anneal', 'NT', 'NX', 'NY', 'W', 'Dp', 'beta');
fprintf('Saved results to: %s\n', results_file);

%% -------------------- FIGURE 1: 2D CONCENTRATION FIELD -------------------------
fig1 = figure('Color', 'w', 'Position', [100 100 1000 850]);
% Convert axes to micrometers for plotting readability.
x_um = x * 1e4;
y_um = y * 1e4;
contourf(x_um, y_um, log10(C), 40, 'LineStyle', 'none');
set(gca, 'YDir', 'reverse');   % y=0 (surface) at top, increasing depth downward
colormap(jet);
cb = colorbar;
cb.Label.String = 'log_{10}( C [cm^{-3}] )';
cb.Label.FontSize = 20;
cb.FontSize = 18;
xlabel('x (\mum)  [across fin]', 'FontSize', 22, 'FontWeight','bold');
ylabel('y (\mum)  [depth into substrate]', 'FontSize', 22, 'FontWeight','bold');
title(sprintf('%s concentration field at t = %.0f s, T = %.0f C', ...
    dopant, t_anneal, T_anneal-273.15), 'FontSize', 20);
set(gca, 'FontSize', 20, 'LineWidth', 1.2);
fig1_path = fullfile(out_dir, 'fem2d_concentration_field.png');
print(fig1, fig1_path, '-dpng', '-r200');
fprintf('Saved figure: %s\n', fig1_path);

%% -------------------- FIGURE 2: CENTERLINE vs 1D REFERENCE ----------------------
fig2 = figure('Color', 'w', 'Position', [100 100 1100 800]);
semilogy(y_um, C_centerline, 'b-', 'LineWidth', 2.5); hold on;
semilogy(y_um, C_1d_ref, 'r--', 'LineWidth', 2.5);
yline(CB, 'k:', 'LineWidth', 1.5);
xlabel('y (\mum)  [depth into substrate]', 'FontSize', 22, 'FontWeight','bold');
ylabel('C (cm^{-3})', 'FontSize', 22, 'FontWeight','bold');
title(sprintf('%s: 2D centerline (x=W/2) vs 1D constant-Ge reference\nmax rel. diff = %.2f%%', ...
    dopant, max_rel_diff_pct), 'FontSize', 20);
legend('2D ADI solution (centerline, laterally-graded + strain-coupled)', ...
       sprintf('1D erfc reference (constant Ge = %.3f, no strain coupling)', x_Ge_avg), ...
       'Background level C_B', 'Location', 'northeast', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 20, 'LineWidth', 1.2);
fig2_path = fullfile(out_dir, 'fem2d_centerline_vs_1d_reference.png');
print(fig2, fig2_path, '-dpng', '-r200');
fprintf('Saved figure: %s\n', fig2_path);

%% -------------------- FIGURE 3: STRAIN FIELD -------------------------------------
fig3 = figure('Color', 'w', 'Position', [100 100 1000 850]);
contourf(x_um, y_um, eps_strain, 40, 'LineStyle', 'none');
set(gca, 'YDir', 'reverse');
colormap(parula);
cb3 = colorbar;
cb3.Label.String = '\epsilon(x,y)  (Vegard-law lattice mismatch strain, dimensionless)';
cb3.Label.FontSize = 18;
cb3.FontSize = 18;
xlabel('x (\mum)  [across fin]', 'FontSize', 22, 'FontWeight','bold');
ylabel('y (\mum)  [depth into substrate]', 'FontSize', 22, 'FontWeight','bold');
title('Lattice-mismatch strain field \epsilon(x,y) from Vegard''s law', 'FontSize', 20);
set(gca, 'FontSize', 20, 'LineWidth', 1.2);
fig3_path = fullfile(out_dir, 'fem2d_strain_field.png');
print(fig3, fig3_path, '-dpng', '-r200');
fprintf('Saved figure: %s\n', fig3_path);

%% -------------------- WALL-CLOCK TIMING ------------------------------------------
elapsed_time = toc(script_start_time);
fprintf('\nTotal wall-clock run time: %.2f s (%.2f min)\n', elapsed_time, elapsed_time/60);
fprintf('=== Done. ===\n');
