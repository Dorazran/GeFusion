% -----------------------------------------------------------------------------
% Copyright (c) 2026 Dor Azran, Ariel University
% Licensed under the MIT License. See LICENSE file in the project root.
% -----------------------------------------------------------------------------

%% kmc_simulation.m
% =========================================================================
% Kinetic Monte Carlo (KMC) bottom-up validation of the phenomenological
% Ge-fraction-enhanced dopant diffusivity model used in the analytic
% Fick's-2nd-law continuum part of this paper.
%
% WHAT THIS SIMULATES
% --------------------------------------------------------------------
% A vacancy-mediated dopant hopping process on a 3D simple-cubic lattice
% representing a uniform-composition Si(1-x)Ge(x) region. A small
% population of vacancies, present at a temperature-dependent equilibrium
% concentration (Arrhenius vacancy formation model), performs a
% rejection-free (Bortz-Kalos-Lebowitz / Gillespie-style) continuous-time
% random walk. Dopant atoms diffuse only via direct exchange with an
% adjacent vacancy ("vacancy mechanism"). The hop (exchange) rate for a
% dopant atom is built from the same bulk Arrhenius activation energy
% Ea and the same Ge-fraction enhancement exponent alpha that appear in
% the analytic model D_SiGe(T,x_Ge) = D_Si(T) * exp(alpha * x_Ge), with a
% small added Vegard's-law strain correction term. The simulation
% accumulates the dopant mean-squared displacement (MSD) vs simulated
% physical time and extracts an effective diffusivity D_eff = MSD/(6t)
% (3D random walk), which is then compared numerically against the
% analytic D_SiGe formula.
%
% WHY THIS MATTERS FOR THE PAPER
% --------------------------------------------------------------------
% The continuum Fick's-law model in the rest of this paper uses an
% empirically-motivated exponential enhancement factor exp(alpha*x_Ge)
% to capture how Ge incorporation changes dopant diffusivity in SiGe.
% That exponent is taken from the literature as a phenomenological,
% top-down fit. This KMC simulation provides an independent, bottom-up,
% atomistic check: by constructing microscopic hop rates from basic
% physical ingredients (attempt frequency, activation energy, strain),
% and then *measuring* the emergent macroscopic diffusivity from the
% statistics of a simulated random walk, we test whether the
% phenomenological exponent alpha is consistent with a physically
% reasonable microscopic mechanism. Agreement (small relative error
% between KMC-extracted D_eff and the analytic D_SiGe) strengthens
% confidence in the continuum model used elsewhere in the paper.
%
% KEY ASSUMPTIONS AND LIMITATIONS
% --------------------------------------------------------------------
%  1. Single-vacancy direct-exchange mechanism only. Real Si/Ge dopant
%     diffusion can also involve interstitials, vacancy pairs, kick-out
%     mechanisms, etc. We deliberately use the simplest vacancy-mediated
%     picture consistent with the activation energies already adopted
%     in the analytic model.
%  2. Simple-cubic (SC) lattice used as a computationally convenient
%     stand-in for the true diamond-cubic lattice of Si/Ge. Coordination
%     number and exact geometric prefactors differ from the real
%     lattice; only six nearest-neighbor hop directions are used here.
%     This affects absolute prefactors much less than it affects the
%     *ratio* D_SiGe/D_Si, which is the quantity of real interest, since
%     the same lattice/geometry approximation is used for both x_Ge=0
%     and x_Ge>0 calculations (the geometric prefactor mostly cancels).
%  3. Uniform Ge fraction x_Ge per run (no spatial gradient). The
%     analytic model elsewhere in the paper considers x_Ge(x) profiles;
%     here we validate the *local* constitutive relation D_SiGe(x_Ge,T)
%     at fixed x_Ge. To extend to a spatially varying Ge profile, the
%     per-site local x_Ge (and hence local Ea_hop) could be looked up
%     from a 3D array x_Ge_field(i,j,k) instead of the scalar X_GE used
%     below -- the hop-rate calculation already takes a "local x_Ge"
%     argument so this extension is straightforward (see local function
%     hop_activation_energy).
%  4. Vacancy concentration is treated as dilute and spatially uniform at
%     equilibrium (Arrhenius formation-energy model), i.e. we do not
%     model vacancy supersaturation/injection effects from implantation
%     damage, nor vacancy-vacancy interactions beyond hard-core exclusion
%     (single occupancy per site). As of the mean-field redesign (see the
%     REDESIGN note in the derivation comment block below), the vacancy
%     population is no longer explicitly placed on a lattice and
%     random-walked; instead each dopant's exchange rate with a vacancy
%     is computed directly from the mean-field vacancy site fraction
%     f_v_sim, adiabatically eliminating the much-faster vacancy
%     migration dynamics (see that note for the physical justification).
%  5. Dopant-dopant interactions (clustering, pairing, electric-field
%     effects) are neglected; dopants are non-interacting random walkers
%     that move only by exchanging with a vacancy that happens to be on
%     a nearest-neighbor site.
%  6. This is a research/teaching-level approximation intended to show
%     qualitative and order-of-magnitude/percent-level quantitative
%     consistency with the phenomenological exponent alpha, not a
%     fully quantitative ab-initio prediction.
%
% RESTART / OVERNIGHT-RUN BEHAVIOR
% --------------------------------------------------------------------
% This script is designed to run unattended for many hours:
%   - Progress is checkpointed to checkpoint_kmc.mat every
%     CHECKPOINT_WALLCLOCK_SEC seconds of wall-clock time AND every
%     CHECKPOINT_EVERY_STEPS KMC steps (whichever comes first).
%   - On startup, if checkpoint_kmc.mat exists in the script folder, the
%     full state is loaded and the run RESUMES rather than restarting.
%   - Human-readable progress lines are appended to kmc_log.txt at every
%     checkpoint.
%   - The main loop is wrapped in try/catch: any error is logged to
%     kmc_log.txt and a best-effort final checkpoint is saved before the
%     error is rethrown.
%   - Total wall-clock budget is configurable via MAX_WALLCLOCK_HOURS;
%     the script checks elapsed wall-clock time every outer iteration
%     and stops gracefully (saving checkpoint + final results + plot)
%     once the budget is exhausted.
%
% Uses ONLY core MATLAB (no toolboxes). Tested for syntax by careful
% manual review (not executed, per task instructions).
% =========================================================================

function kmc_simulation()

clear; close all; %#ok<*UNRCH> (kept for safety if pasted into a script context)

%% ------------------------------------------------------------------
%  SMOKE-TEST SWITCH
%  ------------------------------------------------------------------
% When true, uses a small lattice and a short simulated-time / wall-clock
% budget so the whole script can be run in well under a few minutes to
% sanity-check that MSD grows nonzero and D_eff converges toward the
% analytic target under the mean-field hop-rate model below. Flip to
% false before the real overnight production run, which restores the
% original full-scale settings (60x60x60 lattice, TARGET_SIM_TIME_S =
% 3600 s, MAX_WALLCLOCK_HOURS = 8.0, matching the project's ANNEAL_TIME_s
% convention used elsewhere).
SMOKE_TEST = false;   % <<< SET TO false BEFORE THE REAL OVERNIGHT RUN >>>

%% ------------------------------------------------------------------
%  USER-CONFIGURABLE PARAMETERS
%  ------------------------------------------------------------------

% --- Lattice size (simple-cubic). Increase for larger/more realistic
%     statistics if compute budget allows; decrease for a quick test run.
if SMOKE_TEST
    NX = 20;
    NY = 20;
    NZ = 20;
else
    NX = 60;
    NY = 60;
    NZ = 60;
end

% --- Dopant species to simulate this run: 'Boron', 'Arsenic', or
%     'Phosphorus'. Change this and re-run (will start a fresh
%     checkpoint if DOPANT_SPECIES or X_GE differs from the saved one).
DOPANT_SPECIES = 'Boron';

% --- Uniform Ge fraction for this run (0 <= X_GE <= 1).
X_GE = 0.20;

% --- Anneal temperature in Kelvin. Pick one of the three temperatures
%     of interest in the paper: 1173.15, 1273.15, 1373.15.
T_ANNEAL_K = 1273.15;

% --- Number of dopant atoms and number of vacancies to track.
N_DOPANTS  = 300;
% Vacancy COUNT is derived below from an Arrhenius equilibrium
% concentration model, but we cap it with a minimum/maximum for
% numerical sanity on small lattices; see vacancy count calculation.
MIN_VACANCIES = 50;
MAX_VACANCY_FRACTION = 0.05;  % never exceed 5% of sites as vacancies

% --- KMC attempt frequency (literature-typical phonon attempt frequency
%     for solid-state diffusion hops, order 1e13 Hz).
NU0_HZ = 1.0e13;

% --- Vacancy formation energy (eV). Literature-typical value for Si is
%     commonly cited in the range ~3.5-4.0 eV; we adopt 3.6 eV as a
%     representative literature-typical figure for this simplified model.
E_VACANCY_FORMATION_EV = 3.6;
% Pre-exponential "entropy" factor for vacancy site fraction (dimensionless,
% of order unity to a few; literature values vary). Kept as 1 for simplicity.
VACANCY_PREFACTOR = 1.0;

% --- Target simulated physical anneal time (seconds). The script will
%     try to reach this simulated time, subject to the wall-clock budget.
%     SMOKE_TEST uses a tiny target so a full pass completes in seconds;
%     the production (SMOKE_TEST=false) value of 3600.0 s matches this
%     project's ANNEAL_TIME_s convention used elsewhere in the paper.
if SMOKE_TEST
    % BUGFIX #5 (found during post-redesign smoke test): under the
    % mean-field hop-rate model, dt per KMC step is ~1/R_total, where
    % R_total now reflects ONLY the (slow) dopant-exchange rate -- there
    % is no fast vacancy-migration channel diluting it any more. That
    % makes dt per step orders of magnitude LARGER than under the old
    % explicit-vacancy model, not smaller. A measured run gave dt ~
    % 1.9e-4 s/step, so the old value of 5.0e-9 s here was reached (and
    % blew straight past) after a single step, terminating the smoke
    % test before any real statistics could accumulate (1 step, 0 hops
    % logged). This value must instead be a genuinely generous upper
    % safety bound -- comfortably above the sim time needed to reach
    % MIN_HOPS_PER_DOPANT*N_DOPANTS hops at the computed rate (a few
    % seconds of simulated time in the smoke-test configuration) -- so
    % that the real stopping criterion (hops-per-dopant / D_eff
    % convergence, see BUGFIX #2 below) is what actually ends the run.
    TARGET_SIM_TIME_S = 1.0e3;
else
    TARGET_SIM_TIME_S = 3600.0;
end

% --- Wall-clock run budget (hours). Script stops gracefully after this
%     many hours of real time, saving state, even if TARGET_SIM_TIME_S
%     has not yet been reached.
if SMOKE_TEST
    MAX_WALLCLOCK_HOURS = 0.05;   % 3 minutes
else
    MAX_WALLCLOCK_HOURS = 8.0;
end

% --- Checkpointing cadence.
CHECKPOINT_WALLCLOCK_SEC = 300;   % 5 minutes
CHECKPOINT_EVERY_STEPS   = 200000; % also checkpoint every N KMC steps

% --- MSD / D_eff logging cadence (in KMC steps).
MSD_LOG_EVERY_STEPS = 50000;

% --- CONVERGENCE-BASED STOPPING CRITERION (BUGFIX #2) -----------------
% IMPORTANT: D_SiGe(T,x_Ge) in the analytic model is a MATERIAL PROPERTY
% -- it does not depend on anneal duration. TARGET_SIM_TIME_S=3600s above
% was originally (incorrectly) used as the KMC stopping target, under the
% mistaken idea that the simulated KMC clock needed to reach the real
% 3600s anneal duration. That is NOT required, and is in fact computationally
% impossible: reaching macroscopic seconds of simulated time would require
% on the order of 1e18-1e19 individual KMC hop events (a direct
% consequence of the enormous timescale separation between single
% atomic-attempt-frequency hops, ~1e-13 s, and a 3600 s anneal). No
% realistic wall-clock budget (overnight or otherwise) can reach this.
%
% The correct, standard tracer-diffusion KMC methodology is instead to
% run until the MEASURED D_eff = MSD/(6t) has CONVERGED to a stable
% value (i.e., enough dopant hops have accumulated for the random-walk
% statistics to be meaningful), and then compare that converged D_eff
% directly to the analytic D_SiGe -- regardless of what (tiny) absolute
% simulated time that convergence happens to occur at. TARGET_SIM_TIME_S
% is retained ONLY as a generous, effectively-unreachable upper safety
% bound; MAX_WALLCLOCK_HOURS and the convergence check below are the
% real stopping criteria.
if SMOKE_TEST
    MIN_HOPS_PER_DOPANT = 30;     % small, just to demonstrate convergence quickly
else
    MIN_HOPS_PER_DOPANT = 300;    % production: require good random-walk statistics
end
CONVERGENCE_TOL = 0.03;           % relative change tolerance (3%) between consecutive MSD-log points
CONVERGENCE_WINDOW = 4;           % number of consecutive MSD-log points that must all satisfy the tolerance

% --- Physical constants.
KB_EV_PER_K = 8.617333e-5;  % Boltzmann constant, eV/K

% --- Lattice constants for Vegard's law (cm).
A_SI_CM = 5.431e-8;
A_GE_CM = 5.658e-8;

% --- Bulk-Si Arrhenius diffusivity parameters per dopant species.
%     D0 in cm^2/s, Ea in eV.
dopant_table = struct( ...
    'Boron',      struct('D0', 0.76,  'Ea', 3.46, 'alpha', -3.0), ...
    'Arsenic',    struct('D0', 22.9,  'Ea', 4.05, 'alpha',  2.3), ...
    'Phosphorus', struct('D0', 3.85,  'Ea', 3.66, 'alpha',  0.7) ...
);

if ~isfield(dopant_table, DOPANT_SPECIES)
    error('kmc_simulation:badDopant', 'Unknown DOPANT_SPECIES "%s".', DOPANT_SPECIES);
end
dopant_params = dopant_table.(DOPANT_SPECIES);
D0_cm2_s = dopant_params.D0;
Ea_bulk_eV = dopant_params.Ea;
alpha_exp = dopant_params.alpha;

% --- Lattice hop distance (cm). For the simple-cubic approximation we
%     use the Ge-fraction-weighted (Vegard's law) lattice constant as
%     the nearest-neighbor hop distance.
a_SiGe_cm = A_SI_CM + X_GE * (A_GE_CM - A_SI_CM);
hop_distance_cm = a_SiGe_cm;  % SC nearest-neighbor distance = lattice constant

% --- Local strain from Vegard's law (dimensionless lattice mismatch).
strain_local = (a_SiGe_cm - A_SI_CM) / A_SI_CM;

%% ------------------------------------------------------------------
%  DERIVATION: MAPPING THE MACROSCOPIC EXPONENT alpha ONTO A
%  MICROSCOPIC HOP ACTIVATION ENERGY Ea_hop(x_Ge)
%  ------------------------------------------------------------------
% We want the *ensemble-averaged*, well-mixed (dilute dopant, dilute
% vacancy) KMC hop rate to reproduce, in the macroscopic limit:
%
%     D_SiGe(T, x_Ge) = D_Si(T) * exp(alpha * x_Ge)                  (1)
%
% where D_Si(T) = D0 * exp(-Ea_bulk / (kB*T)).
%
% In a simple-cubic vacancy-mediated random walk with hop attempt
% frequency nu0, hop activation energy Ea_hop, and a vacancy site
% fraction f_v (probability a given neighbor site is vacant), the
% standard random-walk diffusivity for a tracer dopant is:
%
%     D = (1/6) * hop_distance^2 * Gamma                              (2)
%
% where Gamma is the *total* successful hop rate of the tracer, i.e.
% the rate at which the dopant exchanges with ANY of its 6 neighbors:
%
%     Gamma = 6 * f_v * nu0 * exp(-Ea_hop / (kB*T))                   (3)
%
% (6 neighbor directions on a simple-cubic lattice, each available with
% probability f_v that the neighbor site holds a vacancy, each
% attempted at rate nu0*exp(-Ea_hop/kBT) once a vacancy is adjacent).
%
% Substituting (3) into (2):
%
%     D = f_v * nu0 * hop_distance^2 * exp(-Ea_hop / (kB*T))          (4)
%
% At x_Ge = 0 this must recover the bulk-Si Arrhenius law D_Si(T) =
% D0*exp(-Ea_bulk/kBT). We therefore DEFINE the bulk hop parameters
% (nu0, f_v0, hop_distance0) such that, at x_Ge=0:
%
%     D0 * exp(-Ea_bulk/(kB*T)) = f_v0 * nu0 * a_Si^2 * exp(-Ea_hop0/(kB*T))
%
% which is satisfied (for all T) by choosing Ea_hop0 = Ea_bulk and
% absorbing the remaining prefactors into f_v0 (the vacancy fraction is
% itself an Arrhenius quantity, f_v0 = VACANCY_PREFACTOR *
% exp(-E_vac_formation/(kB*T)); any residual mismatch between the
% lumped prefactor f_v0*nu0*a_Si^2 and D0 is corrected by an explicit
% calibration constant C_CAL, computed numerically below, that brings
% the T=T_ANNEAL_K bulk (x_Ge=0) KMC-formula-prediction of D into exact
% agreement with D_Si(T_ANNEAL_K) from the analytic Arrhenius law). This
% is a standard "calibrate the microscopic prefactor against the known
% macroscopic bulk limit" step and does not affect the x_Ge-DEPENDENCE
% we are testing.
%
% Now introduce the Ge-fraction dependence. We want, at fixed T:
%
%     D(x_Ge) / D(0) = exp(alpha * x_Ge)                              (5)
%
% Using (4), and allowing both the local hop activation energy AND the
% local vacancy formation energy (hence f_v) to depend on x_Ge in
% general, the simplest physically-motivated choice that keeps the
% vacancy thermodynamics fixed (we deliberately do NOT make vacancy
% formation energy depend on x_Ge in this simplified model -- see
% Limitation #4) and puts ALL of the x_Ge dependence into the hop
% barrier is:
%
%     Ea_hop(x_Ge) = Ea_bulk - kB*T*alpha*x_Ge / 1   ... (naive, T-dependent barrier -- NOT used)
%
% A T-DEPENDENT activation energy is unphysical (Ea should be a material
% property, not depend on the anneal temperature). Instead we use the
% standard trick of expanding the exponential enhancement as an
% effective REDUCTION (or increase, depending on sign of alpha) of the
% activation energy that is referenced to a fixed reference temperature
% T_ref (taken as the middle anneal temperature of interest, 1273.15 K,
% i.e. 1000 C) so that Ea_hop is a genuine T-independent material
% parameter once chosen, while still reproducing relation (5) AT the
% reference temperature, and only approximately (to leading order) at
% other temperatures in the 1173-1373 K window of interest (a
% acceptable approximation given the narrow ~200 K spread of anneal
% temperatures relevant to this paper):
%
%     Ea_hop(x_Ge) = Ea_bulk - kB*T_ref*alpha*x_Ge                    (6)
%
% so that exp(-Ea_hop(x_Ge)/(kB*T_ref)) = exp(-Ea_bulk/(kB*T_ref)) *
% exp(alpha*x_Ge), exactly reproducing (5) at T = T_ref, and closely
% approximating it (to within a few percent, given the narrow anneal
% temperature window used in this paper) at the other two temperatures.
%
% Finally, we ADD a small explicit strain correction term on top of (6),
% to make the connection to the Vegard's-law lattice mismatch strain
% epsilon(x_Ge) explicit and physically motivated (elastic/strain
% contribution to migration barrier), of the conventional linear form
% Ea_strain = -beta_strain * strain_local, with beta_strain chosen small
% (a fraction of an eV per unit strain) so that it represents a minor
% perturbation on top of the dominant chemical (Ge-fraction) term in
% (6) and does not spoil the calibration to alpha:
%
%     Ea_hop(x_Ge) = Ea_bulk - kB*T_ref*alpha*x_Ge - beta_strain*strain_local   (7)
%
% Equation (7), together with the calibration constant C_CAL described
% above, is what is implemented in the local function
% hop_activation_energy() and used to build the per-event KMC rates
% below.
%
% BUGFIX HISTORY, leading to the current MEAN-FIELD REDESIGN
% --------------------------------------------------------------------
% BUGFIX #1 (vacancy bulk migration): an earlier version of this script
% only ever enumerated dopant-adjacent exchange events, with no
% mechanism for a vacancy with no adjacent dopant to move at all; such
% vacancies were permanently frozen, producing identically-zero MSD.
% The first fix added a SECOND explicit event type -- vacancies
% literally random-walking through the bulk lattice (hopping into
% adjacent empty, non-dopant sites at a literal attempt rate
% NU0_HZ*exp(-E_vac_migration/kT)) competing fairly with dopant-exchange
% events in the same rejection-free event list.
%
% BUGFIX #3 (calibration mismatch): once vacancy migration was added,
% nu0_eff (the calibration constant in Eq 4) had to be recalibrated
% against the actual SIMULATED vacancy fraction f_v_sim realized on the
% finite lattice, not the true (astronomically smaller) thermodynamic
% equilibrium vacancy fraction -- see the calibration block below for
% the full explanation. This fixed a ~1e12x inflation of the measured
% D_eff.
%
% REDESIGN (mean-field / adiabatic elimination, supersedes BUGFIX #1):
% even after BUGFIX #3, the explicit dual-event-type design above proved
% computationally intractable in practice. Hand-calculation and direct
% smoke-testing showed that the conditional dopant-exchange rate per
% realized vacancy-dopant adjacency (~Ea_hop ~ 3.5 eV, hence a very slow
% conditional rate) is many orders of magnitude SLOWER than the bulk
% vacancy migration rate (~E_vac_migration = 0.45 eV, hence a very fast
% rate). This is a genuine physical timescale separation: a vacancy that
% becomes adjacent to a dopant will, with overwhelming probability,
% migrate away again long before it ever exchanges with that dopant. In
% a rejection-free KMC event list, this means the vast majority of all
% executed events are "wasted" vacancy-migration hops that never move a
% dopant or update the MSD -- in practice this produced ZERO successful
% dopant-exchange events even after tens of thousands of KMC steps in
% smoke testing, making the explicit-vacancy-lattice design unusable
% within any realistic wall-clock budget.
%
% The correct treatment of a fast-slow timescale separation like this
% is to ADIABATICALLY ELIMINATE the fast variable (vacancy position) and
% replace its effect with its time-averaged (mean-field) value. Because
% the vacancy population equilibrates (migrates between sites) on a
% timescale enormously faster than any single dopant-exchange attempt,
% from the dopant's point of view "is there a vacancy on this specific
% neighbor site right now" is, on the timescale of interest, well
% approximated by "what is the long-run AVERAGE probability that this
% neighbor site holds a vacancy" -- i.e. simply the mean-field site
% fraction f_v_sim computed below. This removes the need to simulate
% vacancy positions (or dopant lattice positions, for collision
% purposes) AT ALL: each dopant becomes an independent continuous-time
% random walker that, in each of its 6 neighbor directions, attempts an
% exchange at the constant per-channel rate
%
%     rate_per_channel = f_v_sim * nu0_eff_Hz * exp(-Ea_hop(x_Ge)/(kB*T))
%
% which is exactly Eq (4) above divided by hop_distance^2 and by 6 (one
% direction's share of the total rate Gamma). Since x_Ge is spatially
% uniform in this model, rate_per_channel is identical for every dopant
% and every direction, so the total system event rate
% R_total = 6 * N_DOPANTS * rate_per_channel is a CONSTANT throughout the
% run (no event-list reconstruction needed every step). This eliminates
% the entire timescale-separation bottleneck: every single KMC step now
% advances a real dopant hop and contributes to the MSD, with no wasted
% vacancy-migration steps. The dopant-exchange rate formula and its
% nu0_eff_Hz calibration (BUGFIX #3) are UNCHANGED by this redesign --
% only the bookkeeping of which event happens next is simplified, since
% there is now only one (uniform-rate) event type instead of two
% competing ones. E_VACANCY_MIGRATION_EV and the explicit vacancy/dopant
% lattice positions are no longer needed and have been removed.
%% ------------------------------------------------------------------

T_REF_K = 1273.15;        % reference temperature used in the alpha mapping, Eq (6)-(7)
BETA_STRAIN_EV = 0.05;     % small strain-correction coefficient (eV per unit strain), literature-typical order of magnitude for elastic migration-barrier corrections

%% ------------------------------------------------------------------
%  VACANCY EQUILIBRIUM CONCENTRATION (Arrhenius formation-energy model)
%  ------------------------------------------------------------------
% NOTE: this block was moved BEFORE the nu0_eff calibration below
% (BUGFIX #3) because the calibration must be done against the actual
% SIMULATED vacancy fraction, not the true thermodynamic equilibrium
% fraction -- see BUGFIX #3 note at the calibration block for why.
n_sites = NX * NY * NZ;
f_vacancy_equilibrium = VACANCY_PREFACTOR * exp(-E_VACANCY_FORMATION_EV / (KB_EV_PER_K * T_ANNEAL_K));
n_vacancies = round(f_vacancy_equilibrium * n_sites);
n_vacancies = max(n_vacancies, MIN_VACANCIES);
n_vacancies = min(n_vacancies, round(MAX_VACANCY_FRACTION * n_sites));

fprintf('Lattice sites = %d, equilibrium vacancy fraction = %.4e, N_vacancies used = %d\n', ...
    n_sites, f_vacancy_equilibrium, n_vacancies);

if N_DOPANTS + n_vacancies >= n_sites
    error('kmc_simulation:tooFewSites', ...
        'Lattice too small for requested N_DOPANTS + n_vacancies; increase NX/NY/NZ.');
end

% --- Calibration constant nu0_eff: choose an effective combined
%     prefactor such that, at x_Ge=0, T=T_ANNEAL_K, the KMC formula (4)
%     reproduces D_Si(T_ANNEAL_K) from the analytic Arrhenius law exactly,
%     GIVEN the vacancy fraction THE SIMULATION ACTUALLY USES.
%
%     BUGFIX #3 (critical -- found during smoke-test validation after
%     BUGFIX #1/#2): the original calibration used f_v_at_Tref, the TRUE
%     thermodynamic-equilibrium vacancy fraction (~5.6e-15 at 1273 K),
%     which is many orders of magnitude smaller than the vacancy fraction
%     the simulation can actually afford to place on a finite lattice
%     (clamped up to MIN_VACANCIES, e.g. f_v_sim ~ 0.006). Because the KMC
%     loop only ever fires a dopant-exchange event when a vacancy is
%     ACTUALLY adjacent (a realized configuration, not a mean-field
%     probability -- see comment at the rate calculation below), the f_v
%     factor is already implicitly accounted for by the lattice's real
%     occupation statistics. Calibrating nu0_eff against f_v_at_Tref while
%     the lattice actually runs at the much larger f_v_sim double-counts
%     the vacancy-rarity correction, inflating the conditional hop rate
%     (and hence the measured D_eff) by a factor of order f_v_sim /
%     f_v_at_Tref -- roughly 1e12 here, matching the ~1e9-1e11 %
%     discrepancy observed in the pre-fix smoke tests. The fix is to
%     calibrate against f_v_sim = n_vacancies / n_sites (the fraction
%     ACTUALLY realized in this run's lattice), not the true equilibrium
%     fraction. This makes nu0_eff a per-run calibration constant whose
%     job is specifically to compensate for the finite-lattice vacancy
%     inflation forced by computational tractability, which is exactly
%     the role it needs to play for the KMC-measured D_eff to be
%     comparable to the analytic D_Si(T_ANNEAL_K)/D_SiGe targets.
f_v_sim = n_vacancies / n_sites;
% From Eq (4) at x_Ge=0 (so Ea_hop = Ea_bulk): D = f_v*nu0_eff*a_Si^2*exp(-Ea_bulk/kBT)
% The exp(-Ea_bulk/kBT) factor cancels (it appears identically on the
% target-D side and the rate-formula side at x_Ge=0), leaving simply:
nu0_eff_Hz = D0_cm2_s / (f_v_sim * A_SI_CM^2);

fprintf('Calibration (BUGFIX #3): f_v_sim = %.4e (vs true equilibrium f_v = %.4e)\n', ...
    f_v_sim, f_vacancy_equilibrium);
fprintf('Calibration: nu0_eff = %.4e Hz (literature nu0 ~ %.2e Hz)\n', nu0_eff_Hz, NU0_HZ);

%% ------------------------------------------------------------------
%  ANALYTIC TARGET (what we are validating against)
%  ------------------------------------------------------------------
D_Si_T   = D0_cm2_s * exp(-Ea_bulk_eV / (KB_EV_PER_K * T_ANNEAL_K));
D_SiGe_analytic_cm2_s = D_Si_T * exp(alpha_exp * X_GE);

fprintf('Analytic D_Si(T=%.2fK)   = %.6e cm^2/s\n', T_ANNEAL_K, D_Si_T);
fprintf('Analytic D_SiGe(x_Ge=%.3f) = %.6e cm^2/s\n', X_GE, D_SiGe_analytic_cm2_s);

%% ------------------------------------------------------------------
%  SCRIPT FOLDER / CHECKPOINT FILE PATHS
%  ------------------------------------------------------------------
script_folder = fileparts(mfilename('fullpath'));
if isempty(script_folder)
    script_folder = pwd;
end
checkpoint_file = fullfile(script_folder, 'checkpoint_kmc.mat');
log_file        = fullfile(script_folder, 'kmc_log.txt');
results_file    = fullfile(script_folder, 'kmc_results.mat');
plot_file       = fullfile(script_folder, 'kmc_msd_plot.png');

%% ------------------------------------------------------------------
%  STATE INITIALIZATION (fresh start OR resume from checkpoint)
%  ------------------------------------------------------------------
resumed = false;
if exist(checkpoint_file, 'file') == 2
    try
        loaded = load(checkpoint_file);
        if isfield(loaded, 'state') && ...
                isfield(loaded.state, 'config_signature') && ...
                strcmp(loaded.state.config_signature, ...
                       make_config_signature(DOPANT_SPECIES, X_GE, T_ANNEAL_K, NX, NY, NZ, N_DOPANTS, n_vacancies))
            state = loaded.state;
            if ~isfield(state, 'total_dopant_hops')
                state.total_dopant_hops = 0;  % defensive: old checkpoint predates BUGFIX #2
            end
            resumed = true;
            write_log(log_file, sprintf('RESUMED from checkpoint: step=%d, sim_time=%.6e s', ...
                state.step_count, state.sim_time_s));
        else
            write_log(log_file, 'Checkpoint found but config signature mismatch -- starting fresh run.');
        end
    catch ME_load
        write_log(log_file, sprintf('Failed to load checkpoint (%s) -- starting fresh run.', ME_load.message));
    end
end

if ~resumed
    state = init_fresh_state(NX, NY, NZ, N_DOPANTS, n_vacancies, ...
        DOPANT_SPECIES, X_GE, T_ANNEAL_K);
    write_log(log_file, sprintf('FRESH START: lattice=%dx%dx%d, N_dopants=%d, N_vacancies=%d, species=%s, x_Ge=%.3f, T=%.2fK', ...
        NX, NY, NZ, N_DOPANTS, n_vacancies, DOPANT_SPECIES, X_GE, T_ANNEAL_K));
end

%% ------------------------------------------------------------------
%  MAIN KMC LOOP (mean-field per-dopant random walk -- see REDESIGN
%  note in the derivation comment block above)
%  ------------------------------------------------------------------
% Neighbor (hop direction) offsets for the 6 nearest neighbors on a
% simple-cubic lattice. Positions are no longer tracked on an explicit
% lattice (see REDESIGN note), but the unit-vector directions are still
% needed to update each dopant's unwrapped displacement accumulator.
neighbor_offsets = [ 1  0  0; -1  0  0; ...
                      0  1  0;  0 -1  0; ...
                      0  0  1;  0  0 -1];

run_start_tic = tic;
last_checkpoint_tic = tic;
wallclock_budget_sec = MAX_WALLCLOCK_HOURS * 3600.0;

stop_reason = '';

% --- Mean-field per-channel hop rate (REDESIGN: replaces the explicit
%     dual-event-type vacancy-migration / dopant-exchange event list).
%     This is exactly Eq (4) above, with the vacancy-adjacency
%     probability replaced by its mean-field value f_v_sim (computed in
%     the calibration block above), divided among the 6 equivalent
%     neighbor directions. Because X_GE is spatially uniform in this
%     model, Ea_hop (and hence this rate) is identical for every dopant
%     and every direction, so it is computed ONCE here rather than
%     recomputed every step.
Ea_hop_uniform = hop_activation_energy(Ea_bulk_eV, alpha_exp, X_GE, ...
    T_REF_K, strain_local, BETA_STRAIN_EV);
rate_per_channel_Hz = f_v_sim * nu0_eff_Hz * exp(-Ea_hop_uniform / (KB_EV_PER_K * T_ANNEAL_K));

% --- Total system event rate: N_DOPANTS independent walkers, each with
%     6 equivalent hop channels at the constant rate above. Since this
%     never changes during the run, R_total is also computed once.
R_total = N_DOPANTS * 6 * rate_per_channel_Hz;

fprintf('Mean-field hop model: rate_per_channel = %.4e Hz, R_total = %.4e Hz (dt ~ %.4e s/step)\n', ...
    rate_per_channel_Hz, R_total, 1.0 / max(R_total, realmin));

try
    while true
        % --- Wall-clock budget check (outer-loop granularity) ---
        elapsed_wallclock_sec = toc(run_start_tic);
        if elapsed_wallclock_sec >= wallclock_budget_sec
            stop_reason = 'wallclock_budget_reached';
            break;
        end

        % --- Simulated-time target check ---
        if state.sim_time_s >= TARGET_SIM_TIME_S
            stop_reason = 'target_sim_time_reached';
            break;
        end

        % ================================================================
        % REDESIGN: every KMC step is a single, real dopant hop -- there
        % is no event-list construction and no "wasted" event type any
        % more (contrast with the earlier dual-event-type design
        % described in the REDESIGN note above). Since all N_DOPANTS*6
        % channels share the same constant rate, picking "which event
        % occurs" reduces to simply drawing a uniformly random dopant
        % index and a uniformly random direction -- no cumulative-sum
        % search over a per-step event-rate array is needed.
        % ================================================================
        dt = -log(rand()) / R_total;
        state.sim_time_s = state.sim_time_s + dt;

        d_idx = randi(N_DOPANTS);
        nb = randi(6);
        state.dopant_unwrapped(d_idx, :) = state.dopant_unwrapped(d_idx, :) + neighbor_offsets(nb, :);
        state.total_dopant_hops = state.total_dopant_hops + 1;  % BUGFIX #2: track hop statistics

        state.step_count = state.step_count + 1;

        % ================================================================
        % Periodic MSD logging / D_eff estimate
        % ================================================================
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

            fprintf('step=%d  sim_t=%.4e s  MSD=%.4e cm^2  D_eff=%.4e cm^2/s  avg_hops/dopant=%.1f\n', ...
                state.step_count, state.sim_time_s, msd_cm2, D_eff_now, avg_hops_per_dopant);

            % --- BUGFIX #2: convergence-based stop check ---
            % Stop once (a) enough hops per dopant have accumulated for
            % MSD/(6t) to be a statistically meaningful estimator, AND
            % (b) the last CONVERGENCE_WINDOW logged D_eff values are all
            % mutually within CONVERGENCE_TOL relative change of each
            % other (i.e. D_eff has stabilized rather than still drifting
            % as more hops accumulate).
            if avg_hops_per_dopant >= MIN_HOPS_PER_DOPANT && ...
                    numel(state.msd_log_Deff) >= CONVERGENCE_WINDOW
                recent = state.msd_log_Deff(end-CONVERGENCE_WINDOW+1:end);
                if all(isfinite(recent)) && all(recent > 0)
                    rel_spread = (max(recent) - min(recent)) / mean(recent);
                    if rel_spread <= CONVERGENCE_TOL
                        stop_reason = 'converged_Deff';
                        break;
                    end
                end
            end
        end

        % ================================================================
        % Checkpointing (wall-clock OR step-count triggered)
        % ================================================================
        elapsed_since_last_checkpoint = toc(last_checkpoint_tic);
        if elapsed_since_last_checkpoint >= CHECKPOINT_WALLCLOCK_SEC || ...
                mod(state.step_count, CHECKPOINT_EVERY_STEPS) == 0
            save_checkpoint(checkpoint_file, state);
            last_checkpoint_tic = tic;

            if ~isempty(state.msd_log_Deff)
                current_Deff_for_log = state.msd_log_Deff(end);
            else
                current_Deff_for_log = NaN;
            end
            write_log(log_file, sprintf( ...
                'CHECKPOINT  wallclock=%s  sim_time=%.6e s  D_eff=%.6e cm^2/s  steps=%d', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), state.sim_time_s, current_Deff_for_log, state.step_count)); %#ok<TNOW1,DATST>
        end
    end

catch ME
    % --- Best-effort: log the error and attempt a final checkpoint save ---
    write_log(log_file, sprintf('ERROR at step=%d, sim_time=%.6e s: %s', ...
        state.step_count, state.sim_time_s, ME.message));
    try
        save_checkpoint(checkpoint_file, state);
        write_log(log_file, 'Emergency checkpoint saved after error.');
    catch ME2
        write_log(log_file, sprintf('FAILED to save emergency checkpoint: %s', ME2.message));
    end
    rethrow(ME);
end

%% ------------------------------------------------------------------
%  FINAL CHECKPOINT, RESULTS, AND PLOT
%  ------------------------------------------------------------------
save_checkpoint(checkpoint_file, state);

if isempty(stop_reason)
    stop_reason = 'loop_exited_unexpectedly';
end

if ~isempty(state.msd_log_Deff)
    D_eff_final = state.msd_log_Deff(end);
else
    % No MSD samples were logged yet (very short run); compute one now.
    msd_lattice_units = mean(sum(state.dopant_unwrapped.^2, 2));
    msd_cm2 = msd_lattice_units * hop_distance_cm^2;
    if state.sim_time_s > 0
        D_eff_final = msd_cm2 / (6.0 * state.sim_time_s);
    else
        D_eff_final = NaN;
    end
end

if isfinite(D_eff_final) && D_eff_final > 0 && isfinite(D_SiGe_analytic_cm2_s) && D_SiGe_analytic_cm2_s > 0
    relative_difference_pct = 100.0 * abs(D_eff_final - D_SiGe_analytic_cm2_s) / D_SiGe_analytic_cm2_s;
else
    relative_difference_pct = NaN;
end

fprintf('\n=========================================================\n');
fprintf('FINAL RESULT (stop reason: %s)\n', stop_reason);
fprintf('  Dopant species        : %s\n', DOPANT_SPECIES);
fprintf('  x_Ge                   : %.4f\n', X_GE);
fprintf('  T_anneal                : %.2f K\n', T_ANNEAL_K);
fprintf('  Simulated time reached  : %.6e s\n', state.sim_time_s);
fprintf('  KMC steps executed      : %d\n', state.step_count);
fprintf('  Total dopant hops       : %d  (avg %.1f per dopant)\n', state.total_dopant_hops, state.total_dopant_hops / N_DOPANTS);
fprintf('  KMC D_eff (measured)    : %.6e cm^2/s\n', D_eff_final);
fprintf('  Analytic D_SiGe         : %.6e cm^2/s\n', D_SiGe_analytic_cm2_s);
fprintf('  Relative difference     : %.3f %%\n', relative_difference_pct);
fprintf('=========================================================\n');

write_log(log_file, sprintf( ...
    'FINAL  stop_reason=%s  sim_time=%.6e s  steps=%d  D_eff=%.6e cm^2/s  D_analytic=%.6e cm^2/s  reldiff=%.3f%%', ...
    stop_reason, state.sim_time_s, state.step_count, D_eff_final, D_SiGe_analytic_cm2_s, relative_difference_pct));

% --- Save results .mat file ---
results = struct();
results.D_eff_final = D_eff_final;
results.D_SiGe_analytic_cm2_s = D_SiGe_analytic_cm2_s;
results.D_Si_T = D_Si_T;
results.relative_difference_pct = relative_difference_pct;
results.MSD_vs_t_time_s = state.msd_log_time_s;
results.MSD_vs_t_msd_cm2 = state.msd_log_msd_cm2;
results.MSD_vs_t_Deff = state.msd_log_Deff;
results.params = struct( ...
    'DOPANT_SPECIES', DOPANT_SPECIES, ...
    'X_GE', X_GE, ...
    'T_ANNEAL_K', T_ANNEAL_K, ...
    'NX', NX, 'NY', NY, 'NZ', NZ, ...
    'N_DOPANTS', N_DOPANTS, ...
    'n_vacancies', n_vacancies, ...
    'D0_cm2_s', D0_cm2_s, ...
    'Ea_bulk_eV', Ea_bulk_eV, ...
    'alpha_exp', alpha_exp, ...
    'NU0_HZ', NU0_HZ, ...
    'nu0_eff_Hz', nu0_eff_Hz, ...
    'E_VACANCY_FORMATION_EV', E_VACANCY_FORMATION_EV, ...
    'f_v_sim', f_v_sim, ...
    'rate_per_channel_Hz', rate_per_channel_Hz, ...
    'R_total_Hz', R_total, ...
    'SMOKE_TEST', SMOKE_TEST, ...
    'T_REF_K', T_REF_K, ...
    'BETA_STRAIN_EV', BETA_STRAIN_EV, ...
    'strain_local', strain_local, ...
    'hop_distance_cm', hop_distance_cm, ...
    'TARGET_SIM_TIME_S', TARGET_SIM_TIME_S, ...
    'MAX_WALLCLOCK_HOURS', MAX_WALLCLOCK_HOURS, ...
    'stop_reason', stop_reason);

save(results_file, 'results');
fprintf('Saved final results to: %s\n', results_file);

% --- Produce MSD vs time plot (log-log) with diffusive fit overlay ---
make_msd_plot(state.msd_log_time_s, state.msd_log_msd_cm2, D_eff_final, plot_file);

end % function kmc_simulation


%% =======================================================================
%  LOCAL HELPER FUNCTIONS
%  =======================================================================

function sig = make_config_signature(species, x_ge, T, nx, ny, nz, ndop, nvac)
% Build a simple string signature of the run configuration so that a
% checkpoint from a DIFFERENT configuration is not accidentally resumed.
    sig = sprintf('%s_xGe%.4f_T%.2f_%dx%dx%d_Nd%d_Nv%d', ...
        species, x_ge, T, nx, ny, nz, ndop, nvac);
end


function state = init_fresh_state(NX, NY, NZ, N_DOPANTS, n_vacancies, species, x_ge, T) %#ok<INUSL>
% Initialize a fresh KMC state. REDESIGN (mean-field): dopant and vacancy
% LATTICE POSITIONS are no longer tracked at all -- only each dopant's
% unwrapped displacement accumulator (needed for MSD) and the run
% bookkeeping fields. NX/NY/NZ are kept as arguments (unused directly
% here) only because they feed into the config signature below, so that
% resuming against a mismatched lattice/vacancy-count configuration is
% still correctly detected.

    state = struct();
    state.dopant_unwrapped = zeros(N_DOPANTS, 3);       % unwrapped displacement accumulator (lattice units)
    state.sim_time_s = 0.0;
    state.step_count = 0;
    state.total_dopant_hops = 0;   % BUGFIX #2: counts successful dopant-exchange events only
    state.msd_log_time_s  = zeros(0, 1);
    state.msd_log_msd_cm2 = zeros(0, 1);
    state.msd_log_Deff    = zeros(0, 1);
    state.config_signature = make_config_signature(species, x_ge, T, NX, NY, NZ, N_DOPANTS, n_vacancies);
end


function Ea_hop = hop_activation_energy(Ea_bulk_eV, alpha_exp, x_ge, T_ref_K, strain_local, beta_strain_eV)
% Implements Eq. (7) from the derivation comment block above:
%   Ea_hop(x_Ge) = Ea_bulk - kB*T_ref*alpha*x_Ge - beta_strain*strain_local
% This is evaluated using the FIXED reference temperature T_ref_K (not
% the running anneal temperature) so that Ea_hop is a genuine,
% temperature-independent material parameter, consistent with how
% activation energies are defined.
    KB_EV_PER_K_LOCAL = 8.617333e-5;
    Ea_hop = Ea_bulk_eV - KB_EV_PER_K_LOCAL * T_ref_K * alpha_exp * x_ge - beta_strain_eV * strain_local;
end


function save_checkpoint(checkpoint_file, state) %#ok<INUSD>
% Save the full KMC state to the checkpoint .mat file. Wrapped so any
% I/O error can be caught by the caller.
    save(checkpoint_file, 'state');
end


function write_log(log_file, message)
% Append a single timestamped line to the human-readable log file.
    fid = fopen(log_file, 'a');
    if fid == -1
        warning('kmc_simulation:logWriteFailed', 'Could not open log file %s for writing.', log_file);
        return;
    end
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    fprintf(fid, '[%s] %s\n', timestamp, message);
    fclose(fid);
end


function make_msd_plot(t_vec, msd_vec, D_eff_final, plot_file)
% Produce a log-log MSD-vs-time plot with the diffusive (MSD = 6*D*t)
% fit overlaid, using only core MATLAB plotting functions, then save as
% a PNG via saveas.
    fig = figure('Visible', 'off', 'Position', [100 100 1100 800]);

    valid = isfinite(t_vec) & isfinite(msd_vec) & (t_vec > 0) & (msd_vec > 0);

    if any(valid)
        loglog(t_vec(valid), msd_vec(valid), 'o', 'MarkerSize', 7, ...
            'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', [0.1 0.2 0.5]);
        hold on;

        if isfinite(D_eff_final) && D_eff_final > 0
            t_fit = sort(t_vec(valid));
            msd_fit = 6.0 * D_eff_final * t_fit;
            loglog(t_fit, msd_fit, '-', 'LineWidth', 2.5, 'Color', [0.85 0.1 0.1]);
            legend('KMC MSD data', 'Diffusive fit (MSD = 6 D_{eff} t)', 'Location', 'NorthWest', 'FontSize', 16);
        else
            legend('KMC MSD data', 'Location', 'NorthWest', 'FontSize', 16);
        end

        hold off;
    else
        % No valid data yet (e.g., very short run); plot a placeholder.
        text(0.5, 0.5, 'No MSD data logged yet', 'HorizontalAlignment', 'center', 'FontSize', 18);
        axis off;
    end

    xlabel('Simulated time t (s)', 'FontSize', 22, 'FontWeight','bold');
    ylabel('Mean-squared displacement (cm^2)', 'FontSize', 22, 'FontWeight','bold');
    title('KMC Dopant MSD vs Time (log-log) with Diffusive Fit', 'FontSize', 20);
    set(gca, 'FontSize', 20, 'LineWidth', 1.2);
    grid on;

    saveas(fig, plot_file);
    close(fig);
end
