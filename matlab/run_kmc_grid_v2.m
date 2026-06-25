function run_kmc_grid_v2()
% -----------------------------------------------------------------------------
% Copyright (c) 2026 Dor Azran, Ariel University
% Licensed under the MIT License. See LICENSE file in the project root.
% -----------------------------------------------------------------------------

%% run_kmc_grid_v2.m
% =========================================================================
% Deeper, SiGe-composition-focused KMC validation sweep (v2).
%
% Why v2 exists: the first sweep (run_kmc_grid.m, 243 runs) finished in
% only 12.6 minutes -- far faster than intended, because the mean-field
% KMC engine's per-step cost is tiny (~500,000 steps/sec measured) and
% the per-run depth (500 dopants x 600 hops/dopant) was too shallow to
% use a meaningful amount of wall-clock time. The user explicitly asked
% (again) for the run to be as deep, precise, and high-quality as
% possible, with the clear emphasis that the whole point is the SiGe
% (Ge-fraction) composition dependence -- "תתמקד בGeSi" / "תעשה את זה
% כמה שיותר מדויק איכותי ומעמיק". v2 responds on both axes:
%
%   (a) FINER x_Ge RESOLUTION (the primary, paper-relevant axis):
%       17 points from 0 to 0.40 in steps of 0.025 (vs 9 points/0.05 in
%       v1) -- this is the composition axis whose slope gives the
%       phenomenological alpha exponent being validated.
%   (b) MUCH DEEPER per-run statistics:
%       N_DOPANTS:           500   -> 2000   (4x)
%       MIN_HOPS_PER_DOPANT:  600   -> 3000   (5x)
%       CONVERGENCE_TOL:      0.015 -> 0.005  (3x stricter)
%       CONVERGENCE_WINDOW:   6     -> 10
%   (c) MORE REPLICATES per condition (true independent RNG seeds):
%       3 -> 8, giving a real mean +/- standard-error per condition
%       instead of a 3-point spread.
%
% Grid: 3 species x 3 temperatures x 17 x_Ge values x 8 seeds = 1224 runs.
% At the measured ~500k steps/sec and ~6e6 steps/run (2000 dopants x
% 3000 hops/dopant, before any early convergence stop), this is sized to
% run for several hours of genuine wall-clock compute -- not minutes --
% while remaining centered on the SiGe composition dependence.
%
% Results stream incrementally to kmc_grid_v2_results.csv (resumable: a
% combo already present in that CSV is skipped on relaunch). A parallel
% master log (kmc_grid_v2_master_log.txt) records start/end of each run.
% =========================================================================

this_dir = fileparts(mfilename('fullpath'));
csv_file = fullfile(this_dir, 'kmc_grid_v2_results.csv');
master_log = fullfile(this_dir, 'kmc_grid_v2_master_log.txt');
runs_dir = fullfile(this_dir, 'kmc_grid_v2_runs');
if ~exist(runs_dir, 'dir'), mkdir(runs_dir); end

species_list = {'Boron', 'Arsenic', 'Phosphorus'};
T_list = [1173.15, 1273.15, 1373.15];
xge_list = 0:0.025:0.40;            % primary axis: 17 points, fine resolution
seed_list = [101, 202, 303, 404, 505, 606, 707, 808];   % 8 independent replicates

N_DOPANTS = 2000;
MIN_HOPS_PER_DOPANT = 3000;
CONVERGENCE_TOL = 0.005;
CONVERGENCE_WINDOW = 10;
MAX_WALLCLOCK_HOURS = 0.5;          % per-run safety cap (not expected to be hit)
TARGET_SIM_TIME_S = 1.0e4;
MSD_LOG_EVERY_STEPS = 50000;
NX = 60; NY = 60; NZ = 60;

if exist(csv_file, 'file') ~= 2
    fid = fopen(csv_file, 'w');
    fprintf(fid, 'species,x_Ge,T_K,seed,D_eff_final,D_SiGe_analytic_cm2_s,relative_difference_pct,stop_reason,step_count,total_dopant_hops,elapsed_wallclock_sec\n');
    fclose(fid);
    completed = {};
else
    existing = readtable(csv_file);
    completed = strcat(existing.species, '_', string(existing.x_Ge), '_', string(existing.T_K), '_', string(existing.seed));
    completed = cellstr(completed);
end

total_runs = numel(species_list) * numel(T_list) * numel(xge_list) * numel(seed_list);
write_master(master_log, sprintf('=== run_kmc_grid_v2 started: %d species x %d temps x %d x_Ge x %d seeds = %d total runs ===', ...
    numel(species_list), numel(T_list), numel(xge_list), numel(seed_list), total_runs));

grid_start_tic = tic;
run_idx = 0;

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
                    fprintf('[%d/%d] SKIP (already completed): %s x_Ge=%.3f T=%.2fK seed=%d\n', ...
                        run_idx, total_runs, species, x_ge, T_K, seed);
                    continue;
                end

                tag = sprintf('%s_xGe%.3f_T%.0f_seed%d', species, x_ge, T_K, seed);
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
                opts.MSD_LOG_EVERY_STEPS = MSD_LOG_EVERY_STEPS;
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
write_master(master_log, sprintf('=== run_kmc_grid_v2 FINISHED: %d runs attempted, total wallclock = %.2f hours ===', total_runs, total_elapsed_hr));
fprintf('\n=== GRID v2 SWEEP COMPLETE: %d runs, %.2f hours total. Results in %s ===\n', total_runs, total_elapsed_hr, csv_file);

end % function run_kmc_grid_v2


function write_master(log_file, message)
    fid = fopen(log_file, 'a');
    if fid == -1, return; end
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    fprintf(fid, '[%s] %s\n', timestamp, message);
    fclose(fid);
    fprintf('%s\n', message);
end
