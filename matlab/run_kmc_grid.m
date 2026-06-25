function run_kmc_grid()
% -----------------------------------------------------------------------------
% Copyright (c) 2026 Dor Azran, Ariel University
% Licensed under the MIT License. See LICENSE file in the project root.
% -----------------------------------------------------------------------------

%% run_kmc_grid.m
% =========================================================================
% Deep, multi-hour KMC validation sweep, built on top of run_kmc_core.m
% (the refactored, explicitly-RNG-seeded version of kmc_simulation.m's
% engine). Driven by the user's explicit requests:
%   1) "תתמקד בGeSi" -- the sweep's PRIMARY axis is the Ge fraction x_Ge,
%      swept at fine (9-point) resolution from pure Si (x_Ge=0) up to
%      x_Ge=0.40, because the whole point of the paper is the SiGe
%      composition dependence of the diffusivity (the alpha exponent).
%      Species and temperature are swept too (so the validation covers
%      the full parameter space used in the paper), but x_Ge resolution
%      is deliberately the finest of the three.
%   2) "כמה שיותר איכותית ועמוקה" / "לא דקות בודדות" -- each individual
%      run is deeper than the original single-condition run (more
%      dopants, more required hops/dopant, tighter convergence
%      tolerance), AND each condition is run with 3 independent RNG
%      seeds (true replicates, made possible by run_kmc_core's explicit
%      rng(seed) call) so we get a genuine mean +/- spread per condition
%      instead of a single point estimate. Total sweep is sized to run
%      for several hours, not minutes.
%
% GRID:
%   species:      Boron, Arsenic, Phosphorus            (3)
%   T_ANNEAL_K:   1173.15, 1273.15, 1373.15              (3)
%   X_GE:         0, 0.05, 0.10, ..., 0.40               (9)  <-- primary axis
%   replicate seeds: 3 independent seeds per condition   (3)
%   => 3 x 3 x 9 x 3 = 243 independent KMC runs
%
% Per-run rigor (deeper than the original validated single run):
%   N_DOPANTS = 500            (was 300)
%   MIN_HOPS_PER_DOPANT = 600  (was 300)
%   CONVERGENCE_TOL = 0.015    (was 0.03 -- twice as strict)
%   CONVERGENCE_WINDOW = 6     (was 4)
%   MAX_WALLCLOCK_HOURS = 0.1  (6 min hard safety cap per run; the
%                               convergence/hops criterion is expected to
%                               end most runs well before this)
%
% Results are appended incrementally, one row per completed run, to
% kmc_grid_results.csv in this folder -- so progress can be inspected at
% any time during the multi-hour run without waiting for completion.
% A master log (kmc_grid_master_log.txt) records start/end of each run.
% Per-run checkpoint/log files are named uniquely so resumed/parallel
% runs cannot clobber each other's state, and so the WHOLE grid can be
% resumed (skipping already-completed rows) if interrupted.
% =========================================================================

this_dir = fileparts(mfilename('fullpath'));
csv_file = fullfile(this_dir, 'kmc_grid_results.csv');
master_log = fullfile(this_dir, 'kmc_grid_master_log.txt');
runs_dir = fullfile(this_dir, 'kmc_grid_runs');
if ~exist(runs_dir, 'dir'), mkdir(runs_dir); end

species_list = {'Boron', 'Arsenic', 'Phosphorus'};
T_list = [1173.15, 1273.15, 1373.15];
xge_list = 0:0.05:0.40;          % primary axis: 9 points, fine resolution
seed_list = [101, 202, 303];     % 3 independent replicates per condition

N_DOPANTS = 500;
MIN_HOPS_PER_DOPANT = 600;
CONVERGENCE_TOL = 0.015;
CONVERGENCE_WINDOW = 6;
MAX_WALLCLOCK_HOURS = 0.1;
TARGET_SIM_TIME_S = 1.0e3;
NX = 60; NY = 60; NZ = 60;

% --- CSV header (write only if file doesn't already exist) ---
if exist(csv_file, 'file') ~= 2
    fid = fopen(csv_file, 'w');
    fprintf(fid, 'species,x_Ge,T_K,seed,D_eff_final,D_SiGe_analytic_cm2_s,relative_difference_pct,stop_reason,step_count,total_dopant_hops,elapsed_wallclock_sec\n');
    fclose(fid);
    completed = {};
else
    % Resume support: read which (species,x_Ge,T,seed) combos are already done
    existing = readtable(csv_file);
    completed = strcat(existing.species, '_', string(existing.x_Ge), '_', string(existing.T_K), '_', string(existing.seed));
    completed = cellstr(completed);
end

write_master(master_log, sprintf('=== run_kmc_grid started: %d species x %d temps x %d x_Ge x %d seeds = %d total runs ===', ...
    numel(species_list), numel(T_list), numel(xge_list), numel(seed_list), ...
    numel(species_list)*numel(T_list)*numel(xge_list)*numel(seed_list)));

grid_start_tic = tic;
run_idx = 0;
total_runs = numel(species_list) * numel(T_list) * numel(xge_list) * numel(seed_list);

for si = 1:numel(species_list)
    species = species_list{si};
    for ti = 1:numel(T_list)
        T_K = T_list(ti);
        for xi = 1:numel(xge_list)
            x_ge = xge_list(xi);
            for ri = 1:numel(seed_list)
                seed = seed_list(ri);
                run_idx = run_idx + 1;

                key = sprintf('%s_%s_%s_%d', species, num2str(x_ge), num2str(T_K), seed);
                if ismember(key, completed)
                    fprintf('[%d/%d] SKIP (already completed): %s x_Ge=%.2f T=%.2fK seed=%d\n', ...
                        run_idx, total_runs, species, x_ge, T_K, seed);
                    continue;
                end

                tag = sprintf('%s_xGe%.2f_T%.0f_seed%d', species, x_ge, T_K, seed);
                fprintf('\n[%d/%d] RUNNING: %s (elapsed so far: %.1f min)\n', run_idx, total_runs, tag, toc(grid_start_tic)/60);
                write_master(master_log, sprintf('[%d/%d] START %s', run_idx, total_runs, tag));

                opts = struct();
                opts.species = species;
                opts.X_GE = x_ge;
                opts.T_ANNEAL_K = T_K;
                opts.rng_seed = seed;
                opts.N_DOPANTS = N_DOPANTS;
                opts.NX = NX; opts.NY = NY; opts.NZ = NZ;
                opts.MIN_HOPS_PER_DOPANT = MIN_HOPS_PER_DOPANT;
                opts.CONVERGENCE_TOL = CONVERGENCE_TOL;
                opts.CONVERGENCE_WINDOW = CONVERGENCE_WINDOW;
                opts.TARGET_SIM_TIME_S = TARGET_SIM_TIME_S;
                opts.MAX_WALLCLOCK_HOURS = MAX_WALLCLOCK_HOURS;
                opts.checkpoint_file = fullfile(runs_dir, [tag '_checkpoint.mat']);
                opts.log_file = fullfile(runs_dir, [tag '_log.txt']);
                opts.verbose = true;

                try
                    res = run_kmc_core(opts);
                    fid = fopen(csv_file, 'a');
                    fprintf(fid, '%s,%.4f,%.2f,%d,%.6e,%.6e,%.4f,%s,%d,%d,%.2f\n', ...
                        species, x_ge, T_K, seed, res.D_eff_final, res.D_SiGe_analytic_cm2_s, ...
                        res.relative_difference_pct, res.stop_reason, res.step_count, ...
                        res.total_dopant_hops, res.elapsed_wallclock_sec);
                    fclose(fid);
                    write_master(master_log, sprintf('[%d/%d] DONE %s reldiff=%.3f%% stop=%s wallclock=%.1fs', ...
                        run_idx, total_runs, tag, res.relative_difference_pct, res.stop_reason, res.elapsed_wallclock_sec));
                catch ME
                    write_master(master_log, sprintf('[%d/%d] FAILED %s : %s', run_idx, total_runs, tag, ME.message));
                    fid = fopen(csv_file, 'a');
                    fprintf(fid, '%s,%.4f,%.2f,%d,NaN,NaN,NaN,ERROR:%s,0,0,0\n', species, x_ge, T_K, seed, strrep(ME.message, ',', ';'));
                    fclose(fid);
                end
            end
        end
    end
end

total_elapsed_hr = toc(grid_start_tic) / 3600.0;
write_master(master_log, sprintf('=== run_kmc_grid FINISHED: %d runs attempted, total wallclock = %.2f hours ===', total_runs, total_elapsed_hr));
fprintf('\n=== GRID SWEEP COMPLETE: %d runs, %.2f hours total. Results in %s ===\n', total_runs, total_elapsed_hr, csv_file);

end % function run_kmc_grid


function write_master(log_file, message)
    fid = fopen(log_file, 'a');
    if fid == -1, return; end
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    fprintf(fid, '[%s] %s\n', timestamp, message);
    fclose(fid);
    fprintf('%s\n', message);
end
