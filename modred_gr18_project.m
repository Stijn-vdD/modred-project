clearvars;

% Beam specifications
% Double clamped
L   = 5;        % [m] Length of beam
A   = 4e-4;     % [m^2] Area of cross-section
E   = 69e9;     % [N m^-2] Young's Modulus
I   = 1.33e-8;  % [m^4] area moment of inertia
rho = 2700;     % [kg m^-3] Density
L_0 = 4;        % [m] Centre point of force application
W   = 0.05;     % [m] Diameter of force actuator

grid_size = 5000; % Grid size for spatial discretization
t_sim = 0.4;     % [s] Simulation time for time-domain experiments
t_step = 0.001;  % [s] Time step for time-domain experiments

savefigs = true;

%% Question 4: Numerical Solution
% Solve transcendental equation for natural frequencies and mode shapes

N = 12;
x = linspace(0, L, grid_size);
syms beta

% Transcendental equation
eqn_beta = cos(beta*L)*cosh(beta*L) - 1 == 0;

% Solve for spatial eigenvalues
beta_ = zeros(N, 1);
for i = 2:N+1
    beta_(i-1) = double( vpasolve(eqn_beta, beta, i*pi/L) );
end

% Natural frequencies
omega_ = beta_.^2 * sqrt(E*I/(rho*A));

% Mode shape function
w_hat_ = @(x,n) double( ...
    vpa(sin(beta_(n)*x)) - vpa(sinh(beta_(n)*x)) - ...
    ( (vpa(sin(beta_(n)*L)) - vpa(sinh(beta_(n)*L))) / ...
    (vpa(cos(beta_(n)*L)) - vpa(cosh(beta_(n)*L))) ) * ...
    (vpa(cos(beta_(n)*x)) - vpa(cosh(beta_(n)*x))) );

q4_modes = figure(1);clf;
t4 = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:4
    nexttile;
    hold on;
    for j = 1:3
        mode_idx = (i-1)*3 + j;
        if mode_idx <= N
            plot(x, w_hat_(x, mode_idx)*1e3, 'DisplayName', sprintf('Mode %d', mode_idx), 'LineWidth', 1.5);
        end
    end
    hold off;
    grid on;
    xlabel('Position x [m]');
    ylabel('Displacement w(x) [mm]');
    legend('Location', 'southwest', 'Interpreter', 'latex');
end
if savefigs
    uniformFigureStyle(q4_modes, 'Q4_Beam_Modes', 20, 1/2);
end

% Display frequency and natural frequency results
fprintf('Mode\tBeta_n [1/m]\tOmega_n [rad/s]\n');
for n = 1:N
    fprintf('%d\t%.4f\t\t%.4f\n', n, beta_(n), omega_(n));
end

%% Question 5: Initial Deflection Simulation
% Simulate free vibration from initial deflection using modal expansion

% Initial deflection (satisfies clamped BCs)
w0 = @(x) (1 * ((x/L).*(1 - x/L)).^3 .* sin(3*pi*x/L));

t = 0:t_step:t_sim;
N_values = [1, 3, 5];

w_xt = zeros(length(x), length(t),length(N_values));

for idx = 1:length(N_values)
    N = N_values(idx);

    % Project initial deflection onto modal basis
    C_n = zeros(N, 1);
    for n = 1:N
        integrand = @(x) w0(x) .* w_hat_(x, n);
        C_n(n) = (1 / (rho * A * L)) * integral(integrand, 0, L);
    end

    % Time evolution
    for n = 1:N
        w_xt(:, :, idx) = w_xt(:, :, idx) + C_n(n) * w_hat_(x, n)' * cos(omega_(n) * t);
    end
end

% Plot results over time
q5_deflection = figure(2);clf;
for idx = 1:length(N_values)
    N = N_values(idx);

    subplot(length(N_values), 1, idx);
    imagesc(t, x, w_xt(:, :, idx)*1e3);
    axis xy;
    cb = colorbar;
    ylabel(cb, 'Deflection w(x,t) [mm]','Rotation',270,'FontSize',10,'interpreter','latex');
    xlabel('Time [s]');
    ylabel('Position x [m]');
    title(sprintf('N=%d', N));
end

% Compare initial deflection to modal approximation
q5_initial = figure(3);clf;
hold on;
plot(x, w0(x)*1e3, 'k--', 'DisplayName', 'Initial Deflection $w_0(x)$', 'LineWidth', 1.5);
for idx = 1:length(N_values)
    plot(x, w_xt(:,1,idx)*1e3, 'DisplayName', sprintf('Modal Approx. N=%d', ...
        N_values(idx)), 'LineWidth', 1.5);
end
hold off;
grid on;
xlabel('Position x [m]');
ylabel('Deflection $w(x,0)$ [mm]', 'Interpreter', 'latex');
legend('Location', 'southwest', 'Interpreter', 'latex');

if savefigs
    uniformFigureStyle(q5_deflection, 'Q5_Deflection_Simulation', 18, 3/4);
    uniformFigureStyle(q5_initial, 'Q5_Initial_Comparison', 15, 1/2);
end

%% Question 8: Simulation Experiments with State-Space Model
% Build and simulate state-space models using natural mode basis

N_max_q8 = 20;

% Extend beta and omega arrays if needed
if length(beta_) < N_max_q8
    for i = length(beta_)+1:N_max_q8
        beta_(i) = double(vpasolve(eqn_beta, beta, (i+1)*pi/L));
    end
    omega_ = beta_.^2 * sqrt(E*I/(rho*A));
    % Redefine w_hat_ function to use the extended beta_ array
    w_hat_ = @(x,n) double( ...
        vpa(sin(beta_(n)*x)) - vpa(sinh(beta_(n)*x)) - ...
        ( (vpa(sin(beta_(n)*L)) - vpa(sinh(beta_(n)*L))) / ...
        (vpa(cos(beta_(n)*L)) - vpa(cosh(beta_(n)*L))) ) * ...
        (vpa(cos(beta_(n)*x)) - vpa(cosh(beta_(n)*x))) );
end

fprintf('Pre-computing natural mode shapes on grid...\n');
tic;
W_hat_grid = zeros(length(x), N_max_q8);
for n = 1:N_max_q8
    W_hat_grid(:, n) = w_hat_(x, n);
end

% Normalize modes
norm_const_q8 = zeros(N_max_q8, 1);
for n = 1:N_max_q8
    norm_const_q8(n) = sqrt(trapz(x, W_hat_grid(:, n).^2));
end

Phi_nat = W_hat_grid ./ norm_const_q8';
% Function to evaluate mode shape at arbitrary point
phi_q8 = @(x_val,n) interp1(x, Phi_nat(:, n), x_val, 'spline');

% Spatial distribution of input force (actuator location)
ell = @(x_val) double((x_val >= L_0 - W/2) & (x_val <= L_0 + W/2));

% Compute modal input coefficients
ell_grid = ell(x)';
b_n = zeros(N_max_q8, 1);
for n = 1:N_max_q8
    b_n(n) = (1/(rho*A)) * trapz(x, ell_grid .* Phi_nat(:, n));
end

% Precompute mode shapes at beam midpoint for efficiency
x_mid = L/2;
phi_mid = interp1(x, Phi_nat, x_mid, 'spline');

fprintf('  Completed in %.2f sec\n', toc);

% --- Experiment 1: Step Input Response ---
% Test convergence of reduced models with step input
u_step = @(t) 100;
N_values_q8 = [1, 2, 3, 5, 10];

q8_step = figure(4); clf;
hold on;

for idx = 1:length(N_values_q8)
    N_q8 = N_values_q8(idx);
    fprintf('  N = %d modes...', N_q8);
    tic;

    % Build state-space model
    Omega_sq = diag(omega_(1:N_q8).^2);
    A_ss = [zeros(N_q8), eye(N_q8); -Omega_sq, zeros(N_q8)];
    B_ss = [zeros(N_q8, 1); b_n(1:N_q8)];
    x0_ss = zeros(2*N_q8, 1);

    odefun = @(t, x) A_ss*x + B_ss*u_step(t);

    [t_sol, x_sol] = ode45(odefun, t, x0_ss);
    a_coeff = x_sol(:, 1:N_q8);

    % Reconstruct midpoint deflection
    w_mid = a_coeff * phi_mid(1:N_q8)';

    fprintf(' %.2f sec\n', toc);

    if idx == length(N_values_q8)
        plot(t_sol, w_mid*1e3, '--', 'LineWidth', 2, ...
            'DisplayName', sprintf('N=%d', N_q8));
    else
        plot(t_sol, w_mid*1e3, 'LineWidth', 2, ...
            'DisplayName', sprintf('N=%d', N_q8));
    end
end

hold off;
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Deflection at $x=L/2$ [mm]', 'Interpreter', 'latex');
legend('Location', 'southeast', 'Interpreter', 'latex');

if savefigs
    uniformFigureStyle(q8_step, 'Q8_Step_Response', 18, 1/4);
end

% --- Experiment 2: Sinusoidal Input at Resonance ---
% Excite beam at 3rd natural frequency to show resonance
tic;

N_q8 = 10;
omega_resonance = omega_(3);
u_sin = @(t) 50 * sin(omega_resonance * t);

Omega_sq = diag(omega_(1:N_q8).^2);
A_ss = [zeros(N_q8), eye(N_q8); -Omega_sq, zeros(N_q8)];
B_ss = [zeros(N_q8, 1); b_n(1:N_q8)];
x0_ss = zeros(2*N_q8, 1);

odefun = @(t, x) A_ss*x + B_ss*u_sin(t);
[t_sol, x_sol] = ode45(odefun, t, x0_ss);
a_coeff = x_sol(:, 1:N_q8);

t_snapshots = t_sim * [0.25 0.50 0.75 1];

q8_resonance = figure(5); clf;
hold on;

for idx = 1:length(t_snapshots)
    [~, t_idx] = min(abs(t_sol - t_snapshots(idx)));
    w_spatial = Phi_nat(:, 1:N_q8) * a_coeff(t_idx, :)';
    plot(x, w_spatial*1e3, 'LineWidth', 2, ...
        'DisplayName', sprintf('t=%.2fs', t_snapshots(idx)));
end

hold off;
grid on;
xlabel('Position $x$ [m]', 'Interpreter', 'latex');
ylabel('Deflection $w(x,t)$ [mm]', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex');

fprintf('  Completed in %.2f sec\n', toc);

if savefigs
    uniformFigureStyle(q8_resonance, 'Q8_Resonance_Spatial', 18, 1/4);
end

% --- Experiment 3: Free Vibration from Initial Deflection ---
% Compare convergence for free vibration (no input)
w0_q8 = @(x) sin(3*pi*x/L) .* (x/L).^3 .* (1 - x/L).^3;
N_values_q8_2 = [1, 2, 3, 5, 10];

q8_initial = figure(6); clf; hold on;

for idx = 1:length(N_values_q8_2)
    N_q8 = N_values_q8_2(idx);
    fprintf('  N = %d modes...', N_q8);
    tic;

    Omega_sq = diag(omega_(1:N_q8).^2);
    A_ss = [zeros(N_q8), eye(N_q8); -Omega_sq, zeros(N_q8)];
    B_ss = [zeros(N_q8, 1); b_n(1:N_q8)];

    a0_q8 = zeros(N_q8, 1);
    w0_grid = w0_q8(x);
    for n = 1:N_q8
        a0_q8(n) = trapz(x, w0_grid' .* Phi_nat(:, n));
    end

    x0_ss = [a0_q8; zeros(N_q8, 1)];

    odefun = @(t, x) A_ss*x;
    [t_sol, x_sol] = ode45(odefun, t, x0_ss);
    a_coeff = x_sol(:, 1:N_q8);
    w_mid = a_coeff * phi_mid(1:N_q8)';

    fprintf(' %.2f sec\n', toc);

    if idx == length(N_values_q8_2)
        plot(t_sol, w_mid*1e3, '--', 'LineWidth', 2, ...
            'DisplayName', sprintf('N=%d', N_q8));
    else
        plot(t_sol, w_mid*1e3, 'LineWidth', 2, ...
            'DisplayName', sprintf('N=%d', N_q8));
    end
end

grid on;
ylabel('$w(L/2,t)$ [mm]', 'Interpreter', 'latex');
xlabel('Time [s]', 'Interpreter', 'latex');
legend('Location', 'southeast', 'Interpreter', 'latex');

if savefigs
    uniformFigureStyle(q8_initial, 'Q8_Initial_Deflection', 18, 1/4);
end

%% Question 9: POD Basis Computation
% Generate snapshot data and compute POD basis via SVD

% Reference simulation with all modes for snapshot generation
N_ref = N_max_q8;
t_pod = linspace(0, t_sim, 2*grid_size);

Omega_sq_ref = diag(omega_(1:N_ref).^2);
A_ref = [zeros(N_ref), eye(N_ref); -Omega_sq_ref, zeros(N_ref)];
B_ref = [zeros(N_ref, 1); b_n(1:N_ref)];

% Combined input: step + sinusoidal at omega_1 to excite multiple modes
u_combined = @(t) 50 + 30*sin(omega_(1)*t) + 20*sin(omega_(3)*t);

a0_ref = zeros(N_ref, 1);
w0_grid_ref = w0_q8(x);
for n = 1:N_ref
    a0_ref(n) = trapz(x, w0_grid_ref' .* Phi_nat(:, n));
end
x0_ref = [a0_ref; zeros(N_ref, 1)];

odefun_ref = @(t, x) A_ref*x + B_ref*u_combined(t);
[t_ref, x_ref] = ode45(odefun_ref, t_pod, x0_ref);
a_ref = x_ref(:, 1:N_ref);

% Snapshot matrix
W = Phi_nat(:, 1:N_ref) * a_ref';

% SVD: W = U * S * V', POD modes are columns of U
[U_pod_raw, S_pod, V_pod] = svd(W, 'econ');

% Normalize POD modes for L2 orthonormality
dx_pod = x(2) - x(1);
U_pod = U_pod_raw / sqrt(dx_pod);

sigma = diag(S_pod);
energy_content = cumsum(sigma.^2) / sum(sigma.^2);

fprintf('  Singular value energy capture:\n');
R_values_display = [1, 2, 3, 5, 10];
for R = R_values_display
    if R <= length(sigma)
        fprintf('    R=%2d: %.4f%% energy\n', R, energy_content(R)*100);
    end
end

% --- Plot: Singular value spectrum ---
q9_spectrum = figure(7); clf;
subplot(1,2,1);
stem(1:min(10, length(sigma)), sigma(1:min(10, length(sigma))), 'LineWidth', 1.5);
grid on;
xlabel('Mode index $r$', 'Interpreter', 'latex');
ylabel('Singular value $\sigma_r$', 'Interpreter', 'latex');
title('POD Singular Value Spectrum', 'Interpreter', 'latex');

subplot(1,2,2);
plot(1:min(10, length(energy_content)), energy_content(1:min(10, length(energy_content)))*100, 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
hold on;
yline(99, 'k--', 'LineWidth', 1);
hold off;
grid on;
xlabel('POD order $R$', 'Interpreter', 'latex');
ylabel('Cumulative energy [\%]', 'Interpreter', 'latex');
title('Energy Capture vs POD Order', 'Interpreter', 'latex');
legend({'Energy', '99\%'}, 'Location', 'southeast', 'Interpreter', 'latex');

if savefigs
    uniformFigureStyle(q9_spectrum, 'Q9_POD_Spectrum', 18, 1/4);
end

%% Question 10: POD Model Validation
% Build and validate POD-based reduced-order models

R_max = min(10, size(U_pod, 2));

% Create interpolation functions for POD modes
phi_pod = cell(R_max, 1);
for r = 1:R_max
    phi_pod{r} = @(x_val) interp1(x, U_pod(:, r), x_val, 'spline', 0);
end

% Construct second derivative matrix using finite differences
n_pod = length(x);
D2_pod = zeros(n_pod, n_pod);

% 5-point stencil for interior points
for i = 3:n_pod-2
    D2_pod(i, i-2:i+2) = [-1, 16, -30, 16, -1] / (12*dx_pod^2);
end
D2_pod(2, 1:4) = [2, -5, 4, -1] / dx_pod^2;
D2_pod(n_pod-1, n_pod-3:n_pod) = [-1, 4, -5, 2] / dx_pod^2;
D2_pod(1, 1:4) = [2, -5, 4, -1] / dx_pod^2;
D2_pod(n_pod, n_pod-3:n_pod) = [-1, 4, -5, 2] / dx_pod^2;

% Compute stiffness matrix
D2_U = D2_pod * U_pod(:, 1:R_max);
K_pod = (E*I / (rho*A)) * (D2_U' * D2_U * dx_pod);

% Compute POD modal input coefficients
b_pod = zeros(R_max, 1);
for r = 1:R_max
    integrand = zeros(size(x));
    for i = 1:length(x)
        integrand(i) = ell(x(i)) * U_pod(i, r);
    end
    b_pod(r) = (1/(rho*A)) * trapz(x, integrand);
end

% --- Experiment 1: Step Input Response with POD model ---
% Compare POD model convergence with different orders R
fprintf('Experiment 1: Step input response\n');

u_step = @(t) 100;
R_values_pod = [1, 2, 10];

q10_step = figure(8); clf; hold on;

for idx = 1:length(R_values_pod)
    R_pod = R_values_pod(idx);
    fprintf('  R = %d POD modes...', R_pod);
    tic;

    A_pod = [zeros(R_pod), eye(R_pod); -K_pod(1:R_pod, 1:R_pod), zeros(R_pod)];
    B_pod = [zeros(R_pod, 1); b_pod(1:R_pod)];
    x0_pod = zeros(2*R_pod, 1);

    odefun_pod = @(t, x) A_pod*x + B_pod*u_step(t);
    [t_sol_pod, x_sol_pod] = ode45(odefun_pod, t, x0_pod);
    a_pod_coeff = x_sol_pod(:, 1:R_pod);

    w_mid_pod = zeros(length(t_sol_pod), 1);
    for i = 1:length(t_sol_pod)
        for r = 1:R_pod
            w_mid_pod(i) = w_mid_pod(i) + a_pod_coeff(i, r) * phi_pod{r}(x_mid);
        end
    end

    fprintf(' %.2f sec\n', toc);

    if idx == length(R_values_pod)
        plot(t_sol_pod, w_mid_pod*1e3, '--', 'LineWidth', 2, ...
            'DisplayName', sprintf('R=%d', R_pod));
    else
        plot(t_sol_pod, w_mid_pod*1e3, 'LineWidth', 2, ...
            'DisplayName', sprintf('R=%d', R_pod));
    end
end

hold off;
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Deflection at $x=L/2$ [mm]', 'Interpreter', 'latex');
legend('Location', 'southeast', 'Interpreter', 'latex');

if savefigs
    uniformFigureStyle(q10_step, 'Q10_Step_Response_POD', 18, 1/4);
end

% --- Experiment 2: Sinusoidal Input at Resonance with POD model ---
% Test POD model with resonant excitation
fprintf('Experiment 2: Sinusoidal input at resonance\n');
tic;

R_pod = 10;
omega_resonance = omega_(3);
u_sin = @(t) 50 * sin(omega_resonance * t);

A_pod = [zeros(R_pod), eye(R_pod); -K_pod(1:R_pod, 1:R_pod), zeros(R_pod)];
B_pod = [zeros(R_pod, 1); b_pod(1:R_pod)];
x0_pod = zeros(2*R_pod, 1);

odefun_pod = @(t, x) A_pod*x + B_pod*u_sin(t);
[t_sol_pod, x_sol_pod] = ode45(odefun_pod, t, x0_pod);
a_pod_coeff = x_sol_pod(:, 1:R_pod);

t_snapshots = t_sim * [0.25, 0.50, 0.75, 1];

q10_resonance = figure(9); clf; hold on;

for idx = 1:length(t_snapshots)
    [~, t_idx] = min(abs(t_sol_pod - t_snapshots(idx)));

    w_spatial_pod = zeros(size(x));
    for i = 1:length(x)
        for r = 1:R_pod
            w_spatial_pod(i) = w_spatial_pod(i) + a_pod_coeff(t_idx, r) * phi_pod{r}(x(i));
        end
    end

    plot(x, w_spatial_pod*1e3, 'LineWidth', 2, ...
        'DisplayName', sprintf('t=%.2fs', t_snapshots(idx)));
end

hold off;
grid on;
xlabel('Position $x$ [m]', 'Interpreter', 'latex');
ylabel('Deflection $w(x,t)$ [mm]', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex');

fprintf('  Completed in %.2f sec\n', toc);

if savefigs
    uniformFigureStyle(q10_resonance, 'Q10_Resonance_Spatial_POD', 18, 1/4);
end

% --- Experiment 3: Free Vibration from Initial Deflection with POD model ---
% Test POD model with zero input, initial deflection only
fprintf('Experiment 3: Free vibration from initial deflection\n');

w0_q10 = @(x) sin(3*pi*x/L) .* (x/L).^3 .* (1 - x/L).^3;
R_values_pod_2 = [1, 2, 3, 10];

q10_initial = figure(10); clf; hold on;

for idx = 1:length(R_values_pod_2)
    R_pod = R_values_pod_2(idx);
    fprintf('  R = %d POD modes...', R_pod);
    tic;

    A_pod = [zeros(R_pod), eye(R_pod); -K_pod(1:R_pod, 1:R_pod), zeros(R_pod)];
    B_pod = [zeros(R_pod, 1); b_pod(1:R_pod)];

    a0_pod = zeros(R_pod, 1);
    for r = 1:R_pod
        integrand_pod = zeros(size(x));
        for i = 1:length(x)
            integrand_pod(i) = w0_q10(x(i)) * U_pod(i, r);
        end
        a0_pod(r) = trapz(x, integrand_pod);
    end

    x0_pod = [a0_pod; zeros(R_pod, 1)];

    odefun_pod = @(t, x) A_pod*x;
    [t_sol_pod, x_sol_pod] = ode45(odefun_pod, t, x0_pod);
    a_pod_coeff = x_sol_pod(:, 1:R_pod);

    % Compute deflection at midpoint
    w_mid_pod = zeros(length(t_sol_pod), 1);
    for i = 1:length(t_sol_pod)
        for r = 1:R_pod
            w_mid_pod(i) = w_mid_pod(i) + a_pod_coeff(i, r) * phi_pod{r}(x_mid);
        end
    end

    fprintf(' %.2f sec\n', toc);

    if idx == length(R_values_pod_2)
        plot(t_sol_pod, w_mid_pod*1e3, '--', 'LineWidth', 2, ...
            'DisplayName', sprintf('R=%d', R_pod));
    else
        plot(t_sol_pod, w_mid_pod*1e3, 'LineWidth', 2, ...
            'DisplayName', sprintf('R=%d', R_pod));
    end
end

hold off;
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$w(L/2,t)$ [mm]', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex');

if savefigs
    uniformFigureStyle(q10_initial, 'Q10_Initial_Deflection_POD', 18, 1/4);
end

%% --- Comparison: Natural Modes vs POD ---
% Compare accuracy of natural mode basis vs POD basis for same initial condition
fprintf('\nComparison: Natural mode basis vs POD basis\n');

N_compare = 4;
R_compare = 3;
R_truth = 10;  % High-order POD model as reference

Omega_sq_cmp = diag(omega_(1:N_compare).^2);
A_nat = [zeros(N_compare), eye(N_compare); -Omega_sq_cmp, zeros(N_compare)];
B_nat = [zeros(N_compare, 1); b_n(1:N_compare)];

a0_nat = zeros(N_compare, 1);
w0_grid_cmp = w0_q10(x);
for n = 1:N_compare
    a0_nat(n) = trapz(x, w0_grid_cmp' .* Phi_nat(:, n));
end
x0_nat = [a0_nat; zeros(N_compare, 1)];

odefun_nat = @(t, x) A_nat*x;
[t_nat, x_nat] = ode45(odefun_nat, [0, t_sim], x0_nat);

A_pod_cmp = [zeros(R_compare), eye(R_compare); -K_pod(1:R_compare, 1:R_compare), zeros(R_compare)];
a0_pod_cmp = zeros(R_compare, 1);
for r = 1:R_compare
    integrand_pod = zeros(size(x));
    for i = 1:length(x)
        integrand_pod(i) = w0_q10(x(i)) * U_pod(i, r);
    end
    a0_pod_cmp(r) = trapz(x, integrand_pod);
end
x0_pod_cmp = [a0_pod_cmp; zeros(R_compare, 1)];

odefun_pod_cmp = @(t, x) A_pod_cmp*x;
[t_pod_cmp, x_pod_cmp] = ode45(odefun_pod_cmp, [0, t_sim], x0_pod_cmp);

A_pod_tru = [zeros(R_truth), eye(R_truth); -K_pod(1:R_truth, 1:R_truth), zeros(R_truth)];
a0_pod_tru = zeros(R_truth, 1);
for r = 1:R_truth
    integrand_pod = zeros(size(x));
    for i = 1:length(x)
        integrand_pod(i) = w0_q10(x(i)) * U_pod(i, r);
    end
    a0_pod_tru(r) = trapz(x, integrand_pod);
end
x0_pod_tru = [a0_pod_tru; zeros(R_truth, 1)];

odefun_pod_tru = @(t, x) A_pod_tru*x;
[t_pod_tru, x_pod_tru] = ode45(odefun_pod_tru, [0, t_sim], x0_pod_tru);

w_mid_nat = x_nat(:, 1:N_compare) * phi_mid(1:N_compare)';

w_mid_pod_cmp = zeros(length(t_pod_cmp), 1);
for i = 1:length(t_pod_cmp)
    for r = 1:R_compare
        w_mid_pod_cmp(i) = w_mid_pod_cmp(i) + x_pod_cmp(i, r) * phi_pod{r}(x_mid);
    end
end

w_mid_pod_tru = zeros(length(t_pod_tru), 1);
for i = 1:length(t_pod_tru)
    for r = 1:R_truth
        w_mid_pod_tru(i) = w_mid_pod_tru(i) + x_pod_tru(i, r) * phi_pod{r}(x_mid);
    end
end

q10_comparison = figure(11); clf; hold on;
plot(t_nat, w_mid_nat*1e3, 'LineWidth', 2, 'DisplayName', sprintf('Natural modes (N=%d)', N_compare));
plot(t_pod_cmp, w_mid_pod_cmp*1e3, 'LineWidth', 2, 'DisplayName', sprintf('POD (R=%d)', R_compare));
plot(t_pod_tru, w_mid_pod_tru*1e3, '--', 'LineWidth', 2, 'DisplayName', sprintf('POD (R=%d)', R_truth));
hold off;
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('$w(L/2,t)$ [mm]', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Interpreter', 'latex');

if savefigs
    uniformFigureStyle(q10_comparison, 'Q10_Comparison_Natural_vs_POD', 18, 1/4);
end

%% Figure export function
function result = uniformFigureStyle(figureHandle,fileName,width,aspectRatio)
% Author: S. van den Dungen, 2025

% Change the figure size
figureHandle.Units               = 'centimeters';
figureHandle.Position(3)         = width;
figureHandle.Position(4)         = aspectRatio * width;

% Get all axes that are in the figure
allAxesInFigure = findall(figureHandle,'type','axes');

% Set the figure rendered to painters (vector graphics)
% set(figureHandle, 'Renderer', 'painters');

result = 0;

for k = 1:numel(allAxesInFigure)

    ax = allAxesInFigure(k);

    % Set the box line width to 1 pt
    ax.LineWidth = 1;
    % Show the entire box, not just the axes
    ax.Box = 'on';
    % Remove unnnecessay whitespace
    ax.LooseInset = max(ax.TightInset, 0.02);
    axis(ax,'tight')

    % Change the grid
    ax.GridColor = [0.15 0.15 0.15];
    ax.MinorGridColor = [0.1 0.1 0.1];

    ax.FontSize = 10;
    % Set to latex interpreter
    ax.TickLabelInterpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    if isprop(ax,'ZLabel')
        ax.ZLabel.Interpreter = 'latex';
    end
    if isprop(ax,'Title')
        ax.Title.Interpreter = 'latex';
    end

    result = result + 1;
end

% 'Remove' all borders from legends
leg = findall(figureHandle, 'Tag', 'legend');
if numel(leg)
    % Set the edge color to white ,without removing its white
    % background
    set(leg, 'EdgeColor', [1 1 1]);
end

% Make output folders if necessary
if ~(exist('outputPNG','dir'))
    mkdir('./outputPNG')
end
if ~exist('outputEPS','dir')
    mkdir('./outputEPS')
end
% if ~exist('outputFIG','dir')
%     mkdir('./outputFIG')
% end
% Ensure the exported files have the same size as the figure
figureHandle.PaperPositionMode = 'auto';

drawnow
pause(0.5) % Ensure the figure is fully rendered before exporting

% Export to files
exportgraphics(figureHandle,['outputEPS/',fileName,'.eps'], 'ContentType','vector', 'BackgroundColor','none', 'Resolution',600);
exportgraphics(figureHandle,['outputPNG/',fileName,'.png'], 'BackgroundColor','none', 'Resolution',600);
% savefig(figureHandle,['outputFIG/',fileName]);

% results is number of handles we treated.

end