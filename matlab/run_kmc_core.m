function results = run_kmc_core(opts)
% -----------------------------------------------------------------------------
% Copyright (c) 2026 Dor Azran, Ariel University
% Licensed under the MIT License. See LICENSE file in the project root.
% -----------------------------------------------------------------------------

%% run_kmc_core.m
% =========================================================================
% Core mean-field KMC engine, factored out of kmc_simulation.m so that it
% can be called repeatedly with different parameters (species, x_Ge, T,
% replicate RNG seed, ...) by a grid-sweep driver (see run_kmc_grid.m),
% while kmc_simulation.m itself remains a thin single-run wrapper around
% this function for backward-compatible interactive/smoke-test use.
%
% All of the physics, derivation, and bug-history notes documented at
% length in kmc_simulation.m apply unchanged here -- this file is a
% mechanical refactor (parameterization + explicit RNG seeding for
% genuine statistical independence between replicate runs), not a
% physics change. See kmc_simulation.m header comments for the full
% derivation and bugfix history.
%
% INPUT: opts struct with fields (all required unless noted "default"):
%   species              'Boron' | 'Arsenic' | 'Phosphorus'
%   X_GE                 uniform Ge fraction, 0<=X_GE<=1
%   T_ANNEAL_K           anneal temperature, K
%   rng_seed             integer seed for rng() -- REQUIRED for genuine
%                        statistical independence between replicate runs
%                        (kmc_simulation.m's original single-run script
%                        never set this explicitly, so every run used
%                        MATLAB's fresh-process default stream and was
%                        bit-for-bit IDENTICAL to every other run with
%                        the same parameters -- not a meaningful
%                        replicate). Pass a different integer per
%                        replicate to get independent random walks.
%   N_DOPANTS            default 300
%   NX, NY, NZ           default 60,60,60
%   MIN_HOPS_PER_DOPANT  default 300
%   CONVERGENCE_TOL      default 0.03
%   CONVERGENCE_WINDOW   default 4
%   TARGET_SIM_TIME_S    default 1.0e3 (generous upper safety bound -- see
%                        BUGFIX #5 in kmc_simulation.m for why this must
%                        be large under the mean-field rate model)
%   MAX_WALLCLOCK_HOURS  default 0.25 (15 min per run; grid driver caller
%                        should size this so N_total_runs * this stays
%                        within the overall overnight budget)
%   MSD_LOG_EVERY_STEPS  default 50000
%   CHECKPOINT_WALLCLOCK_SEC default 300
%   CHECKPOINT_EVERY_STEPS   default 200000
%   checkpoint_file      path for this run's checkpoint .mat (caller must
%                        give each (species,x_Ge,T,seed) combo its OWN
%                        checkpoint file so resumes can't cross-contaminate)
%   log_file             path for this run's human-readable log
%   verbose              default true; if false, suppresses fprintf (CSV/
%                        log-file output is unaffected)
%
% OUTPUT: results struct (same fields as kmc_simulation.m's results.mat,
% plus elapsed_wallclock_sec and rng_seed for provenance):
%   D_eff_final, D_SiGe_analytic_cm2_s, D_Si_T, relative_difference_pct,
%   MSD_vs_t_time_s, MSD_vs_t_msd_cm2, MSD_vs_t_Deff, params (struct),
%   elapsed_wallclock_sec
% =========================================================================

% --- Fill in defaults for any opts fields not supplied ---
opts = set_default(opts, 'N_DOPANTS', 300);
opts = set_default(opts, 'NX', 60);
opts = set_default(opts, 'NY', 60);
opts = set_default(opts, 'NZ', 60);
opts = set_default(opts, 'MIN_HOPS_PER_DOPANT', 300);
opts = set_default(opts, 'CONVERGENCE_TOL', 0.03);
opts = set_default(opts, 'CONVERGENCE_WINDOW', 4);
opts = set_default(opts, 'TARGET_SIM_TIME_S', 1.0e3);
opts = set_default(opts, 'MAX_WALLCLOCK_HOURS', 0.25);
opts = set_default(opts, 'MSD_LOG_EVERY_STEPS', 50000);
opts = set_default(opts, 'CHECKPOINT_WALLCLOCK_SEC', 300);
opts = set_default(opts, 'CHECKPOINT_EVERY_STEPS', 200000);
opts = set_default(opts, 'verbose', true);

if ~isfield(opts, 'rng_seed')
    error('run_kmc_core:noSeed', ...
        'opts.rng_seed is required (pass a distinct integer per replicate for statistical independence).');
end
rng(opts.rng_seed, 'twister');   % BUGFIX #6: explicit per-run seeding (see header note) -- without
                                  % this, every run (and every "replicate") is bit-for-bit identical.

DOPANT_SPECIES = opts.species;
X_GE = opts.X_GE;
T_ANNEAL_K = opts.T_ANNEAL_K;
N_DOPANTS = opts.N_DOPANTS;
NX = opts.NX; NY = opts.NY; NZ = opts.NZ;
MIN_HOPS_PER_DOPANT = opts.MIN_HOPS_PER_DOPANT;
CONVERGENCE_TOL = opts.CONVERGENCE_TOL;
CONVERGENCE_WINDOW = opts.CONVERGENCE_WINDOW;
TARGET_SIM_TIME_S = opts.TARGET_SIM_TIME_S;
MAX_WALLCLOCK_HOURS = opts.MAX_WALLCLOCK_HOURS;
MSD_LOG_EVERY_STEPS = opts.MSD_LOG_EVERY_STEPS;
CHECKPOINT_WALLCLOCK_SEC = opts.CHECKPOINT_WALLCLOCK_SEC;
CHECKPOINT_EVERY_STEPS = opts.CHECKPOINT_EVERY_STEPS;
checkpoint_file = opts.checkpoint_file;
log_file = opts.log_file;
verbose = opts.verbose;

% --- Fixed physical constants / vacancy model (unchanged from kmc_simulation.m) ---
MIN_VACANCIES = 50;
MAX_VACANCY_FRACTION = 0.05;
NU0_HZ = 1.0e13;
E_VACANCY_FORMATION_EV = 3.6;
VACANCY_PREFACTOR = 1.0;
KB_EV_PER_K = 8.617333e-5;
A_SI_CM = 5.431e-8;
A_GE_CM = 5.658e-8;
T_REF_K = 1273.15;
BETA_STRAIN_EV = 0.05;

dopant_table = struct( ...
    'Boron',      struct('D0', 0.76,  'Ea', 3.46, 'alpha', -3.0), ...
    'Arsenic',    struct('D0', 22.9,  'Ea', 4.05, 'alpha',  2.3), ...
    'Phosphorus', struct('D0', 3.85,  'Ea', 3.66, 'alpha',  0.7) ...
);
if ~isfield(dopant_table, DOPANT_SPECIES)
    error('run_kmc_core:badDopant', 'Unknown species "%s".', DOPANT_SPECIES);
end
dp = dopant_table.(DOPANT_SPECIES);
D0_cm2_s = dp.D0; Ea_bulk_eV = dp.Ea; alpha_exp = dp.alpha;

a_SiGe_cm = A_SI_CM + X_GE * (A_GE_CM - A_SI_CM);
hop_distance_cm = a_SiGe_cm;
strain_local = (a_SiGe_cm - A_SI_CM) / A_SI_CM;

n_sites = NX * NY * NZ;
f_vacancy_equilibrium = VACANCY_PREFACTOR * exp(-E_VACANCY_FORMATION_EV / (KB_EV_PER_K * T_ANNEAL_K));
n_vacancies = round(f_vacancy_equilibrium * n_sites);
n_vacancies = max(n_vacancies, MIN_VACANCIES);
n_vacancies = min(n_vacancies, round(MAX_VACANCY_FRACTION * n_sites));

if N_DOPANTS + n_vacancies >= n_sites
    error('run_kmc_core:tooFewSites', 'Lattice too small for N_DOPANTS + n_vacancies; increase NX/NY/NZ.');
end

f_v_sim = n_vacancies / n_sites;
nu0_eff_Hz = D0_cm2_s / (f_v_sim * A_SI_CM^2);   % BUGFIX #3 calibration, unchanged

D_Si_T = D0_cm2_s * exp(-Ea_bulk_eV / (KB_EV_PER_K * T_ANNEAL_K));
D_SiGe_analytic_cm2_s = D_Si_T * exp(alpha_exp * X_GE);

if verbose
    fprintf('[seed=%d] %s x_Ge=%.3f T=%.2fK : Lattice=%dx%dx%d, f_v_sim=%.4e, nu0_eff=%.4e Hz, D_analytic=%.6e cm^2/s\n', ...
        opts.rng_seed, DOPANT_SPECIES, X_GE, T_ANNEAL_K, NX, NY, NZ, f_v_sim, nu0_eff_Hz, D_SiGe_analytic_cm2_s);
end

% --- Resume-or-fresh state ---
resumed = false;
if exist(checkpoint_file, 'file') == 2
    try
        loaded = load(checkpoint_file);
        sig = make_config_signature(DOPANT_SPECIES, X_GE, T_ANNEAL_K, NX, NY, NZ, N_DOPANTS, n_vacancies, opts.rng_seed);
        if isfield(loaded, 'state') && isfield(loaded.state, 'config_signature') && ...
                strcmp(loaded.state.config_signature, sig)
            state = loaded.state;
            if ~isfield(state, 'total_dopant_hops'), state.total_dopant_hops = 0; end
            resumed = true;
            write_log(log_file, sprintf('RESUMED from checkpoint: step=%d, sim_time=%.6e s', state.step_count, state.sim_time_s));
        else
            write_log(log_file, 'Checkpoint found but config signature mismatch -- starting fresh run.');
        end
    catch ME_load
        write_log(log_file, sprintf('Failed to load checkpoint (%s) -- starting fresh run.', ME_load.message));
    end
end
if ~resumed
    state = init_fresh_state(NX, NY, NZ, N_DOPANTS, n_vacancies, DOPANT_SPECIES, X_GE, T_ANNEAL_K, opts.rng_seed);
    write_log(log_file, sprintf('FRESH START: lattice=%dx%dx%d, N_dopants=%d, N_vacancies=%d, species=%s, x_Ge=%.3f, T=%.2fK, seed=%d', ...
        NX, NY, NZ, N_DOPANTS, n_vacancies, DOPANT_SPECIES, X_GE, T_ANNEAL_K, opts.rng_seed));
end

neighbor_offsets = [ 1  0  0; -1  0  0; 0  1  0;  0 -1  0; 0  0  1;  0  0 -1];

run_start_tic = tic;
last_checkpoint_tic = tic;
wallclock_budget_sec = MAX_WALLCLOCK_HOURS * 3600.0;
stop_reason = '';

Ea_hop_uniform = hop_activation_energy(Ea_bulk_eV, alpha_exp, X_GE, T_REF_K, strain_local, BETA_STRAIN_EV);
rate_per_channel_Hz = f_v_sim * nu0_eff_Hz * exp(-Ea_hop_uniform / (KB_EV_PER_K * T_ANNEAL_K));
R_total = N_DOPANTS * 6 * rate_per_channel_Hz;

try
    while true
        elapsed_wallclock_sec = toc(run_start_tic);
        if elapsed_wallclock_sec >= wallclock_budget_sec
            stop_reason = 'wallclock_budget_reached'; break;
        end
        if state.sim_time_s >= TARGET_SIM_TIME_S
            stop_reason = 'target_sim_time_reached'; break;
        end

        dt = -log(rand()) / R_total;
        state.sim_time_s = state.sim_time_s + dt;
        d_idx = randi(N_DOPANTS);
        nb = randi(6);
        state.dopant_unwrapped(d_idx, :) = state.dopant_unwrapped(d_idx, :) + neighbor_offsets(nb, :);
        state.total_dopant_hops = state.total_dopant_hops + 1;
        state.step_count = state.step_count + 1;

        if mod(state.step_count, MSD_LOG_EVERY_STEPS) == 0
            msd_lattice_units = mean(sum(state.dopant_unwrapped.^2, 2));
            msd_cm2 = msd_lattice_units * hop_distance_cm^2;
            if state.sim_time_s > 0
                D_eff_now = msd_cm2 / (6.0 * state.sim_time_s);
            else
                D_eff_now = NaN;
            end
            state.msd_log_time_s(end+1, 1) = state.sim_time_s;
            state.msd_log_msd_cm2(end+1, 1) = msd_cm2;
            state.msd_log_Deff(end+1, 1) = D_eff_now;
            avg_hops_per_dopant = state.total_dopant_hops / N_DOPANTS;

            if verbose
                fprintf('  [seed=%d] step=%d sim_t=%.4e D_eff=%.4e avg_hops/dopant=%.1f\n', ...
                    opts.rng_seed, state.step_count, state.sim_time_s, D_eff_now, avg_hops_per_dopant);
            end

            if avg_hops_per_dopant >= MIN_HOPS_PER_DOPANT && numel(state.msd_log_Deff) >= CONVERGENCE_WINDOW
                recent = state.msd_log_Deff(end-CONVERGENCE_WINDOW+1:end);
                if all(isfinite(recent)) && all(recent > 0)
                    rel_spread = (max(recent) - min(recent)) / mean(recent);
                    if rel_spread <= CONVERGENCE_TOL
                        stop_reason = 'converged_Deff'; break;
                    end
                end
            end
        end

        elapsed_since_last_checkpoint = toc(last_checkpoint_tic);
        if elapsed_since_last_checkpoint >= CHECKPOINT_WALLCLOCK_SEC || mod(state.step_count, CHECKPOINT_EVERY_STEPS) == 0
            save_checkpoint(checkpoint_file, state);
            last_checkpoint_tic = tic;
            if ~isempty(state.msd_log_Deff), cur = state.msd_log_Deff(end); else, cur = NaN; end
            write_log(log_file, sprintf('CHECKPOINT wallclock=%s sim_time=%.6e D_eff=%.6e steps=%d', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), state.sim_time_s, cur, state.step_count)); %#ok<TNOW1,DATST>
        end
    end
catch ME
    write_log(log_file, sprintf('ERROR at step=%d, sim_time=%.6e s: %s', state.step_count, state.sim_time_s, ME.message));
    try
        save_checkpoint(checkpoint_file, state);
        write_log(log_file, 'Emergency checkpoint saved after error.');
    catch ME2
        write_log(log_file, sprintf('FAILED to save emergency checkpoint: %s', ME2.message));
    end
    rethrow(ME);
end

elapsed_wallclock_sec_final = toc(run_start_tic);
save_checkpoint(checkpoint_file, state);
if isempty(stop_reason), stop_reason = 'loop_exited_unexpectedly'; end

if ~isempty(state.msd_log_Deff)
    D_eff_final = state.msd_log_Deff(end);
else
    msd_lattice_units = mean(sum(state.dopant_unwrapped.^2, 2));
    msd_cm2 = msd_lattice_units * hop_distance_cm^2;
    if state.sim_time_s > 0, D_eff_final = msd_cm2 / (6.0 * state.sim_time_s); else, D_eff_final = NaN; end
end

if isfinite(D_eff_final) && D_eff_final > 0 && isfinite(D_SiGe_analytic_cm2_s) && D_SiGe_analytic_cm2_s > 0
    relative_difference_pct = 100.0 * abs(D_eff_final - D_SiGe_analytic_cm2_s) / D_SiGe_analytic_cm2_s;
else
    relative_difference_pct = NaN;
end

write_log(log_file, sprintf('FINAL stop_reason=%s sim_time=%.6e steps=%d D_eff=%.6e D_analytic=%.6e reldiff=%.3f%% seed=%d', ...
    stop_reason, state.sim_time_s, state.step_count, D_eff_final, D_SiGe_analytic_cm2_s, relative_difference_pct, opts.rng_seed));

if verbose
    fprintf('  -> DONE [%s, x_Ge=%.3f, T=%.2fK, seed=%d]: D_eff=%.6e vs analytic=%.6e (reldiff=%.3f%%, stop=%s, steps=%d, wallclock=%.1fs)\n', ...
        DOPANT_SPECIES, X_GE, T_ANNEAL_K, opts.rng_seed, D_eff_final, D_SiGe_analytic_cm2_s, relative_difference_pct, stop_reason, state.step_count, elapsed_wallclock_sec_final);
end

results = struct();
results.D_eff_final = D_eff_final;
results.D_SiGe_analytic_cm2_s = D_SiGe_analytic_cm2_s;
results.D_Si_T = D_Si_T;
results.relative_difference_pct = relative_difference_pct;
results.MSD_vs_t_time_s = state.msd_log_time_s;
results.MSD_vs_t_msd_cm2 = state.msd_log_msd_cm2;
results.MSD_vs_t_Deff = state.msd_log_Deff;
results.elapsed_wallclock_sec = elapsed_wallclock_sec_final;
results.stop_reason = stop_reason;
results.step_count = state.step_count;
results.total_dopant_hops = state.total_dopant_hops;
results.params = struct( ...
    'DOPANT_SPECIES', DOPANT_SPECIES, 'X_GE', X_GE, 'T_ANNEAL_K', T_ANNEAL_K, ...
    'NX', NX, 'NY', NY, 'NZ', NZ, 'N_DOPANTS', N_DOPANTS, 'n_vacancies', n_vacancies, ...
    'D0_cm2_s', D0_cm2_s, 'Ea_bulk_eV', Ea_bulk_eV, 'alpha_exp', alpha_exp, ...
    'NU0_HZ', NU0_HZ, 'nu0_eff_Hz', nu0_eff_Hz, 'E_VACANCY_FORMATION_EV', E_VACANCY_FORMATION_EV, ...
    'f_v_sim', f_v_sim, 'rate_per_channel_Hz', rate_per_channel_Hz, 'R_total_Hz', R_total, ...
    'T_REF_K', T_REF_K, 'BETA_STRAIN_EV', BETA_STRAIN_EV, 'strain_local', strain_local, ...
    'hop_distance_cm', hop_distance_cm, 'TARGET_SIM_TIME_S', TARGET_SIM_TIME_S, ...
    'MAX_WALLCLOCK_HOURS', MAX_WALLCLOCK_HOURS, 'MIN_HOPS_PER_DOPANT', MIN_HOPS_PER_DOPANT, ...
    'CONVERGENCE_TOL', CONVERGENCE_TOL, 'CONVERGENCE_WINDOW', CONVERGENCE_WINDOW, ...
    'rng_seed', opts.rng_seed, 'stop_reason', stop_reason);

end % function run_kmc_core


function opts = set_default(opts, field, value)
    if ~isfield(opts, field) || isempty(opts.(field))
        opts.(field) = value;
    end
end

function sig = make_config_signature(species, x_ge, T, nx, ny, nz, ndop, nvac, seed)
    sig = sprintf('%s_xGe%.4f_T%.2f_%dx%dx%d_Nd%d_Nv%d_seed%d', species, x_ge, T, nx, ny, nz, ndop, nvac, seed);
end

function state = init_fresh_state(NX, NY, NZ, N_DOPANTS, n_vacancies, species, x_ge, T, seed) %#ok<INUSL>
    state = struct();
    state.dopant_unwrapped = zeros(N_DOPANTS, 3);
    state.sim_time_s = 0.0;
    state.step_count = 0;
    state.total_dopant_hops = 0;
    state.msd_log_time_s  = zeros(0, 1);
    state.msd_log_msd_cm2 = zeros(0, 1);
    state.msd_log_Deff    = zeros(0, 1);
    state.config_signature = make_config_signature(species, x_ge, T, NX, NY, NZ, N_DOPANTS, n_vacancies, seed);
end

function Ea_hop = hop_activation_energy(Ea_bulk_eV, alpha_exp, x_ge, T_ref_K, strain_local, beta_strain_eV)
    KB_EV_PER_K_LOCAL = 8.617333e-5;
    Ea_hop = Ea_bulk_eV - KB_EV_PER_K_LOCAL * T_ref_K * alpha_exp * x_ge - beta_strain_eV * strain_local;
end

function save_checkpoint(checkpoint_file, state) %#ok<INUSD>
    save(checkpoint_file, 'state');
end

function write_log(log_file, message)
    fid = fopen(log_file, 'a');
    if fid == -1
        warning('run_kmc_core:logWriteFailed', 'Could not open log file %s for writing.', log_file);
        return;
    end
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    fprintf(fid, '[%s] %s\n', timestamp, message);
    fclose(fid);
end
