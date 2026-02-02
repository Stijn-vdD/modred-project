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

savefigs = true;

%% Question 4: Numerical Solution
% Numerically solve the transcendental equation for N>=10 or more frequencies
% beta_n, the corresponding natural frequencies omega_n, and the natural
% vibration modes w_hat_n(x).

N = 12;
syms beta

eqn_beta = cos(beta*L)*cosh(beta*L) - 1 == 0;

beta_ = zeros(N, 1);
for i = 2:N+1
    beta_(i-1) = double( vpasolve(eqn_beta, beta, i*pi/L) );
end

omega_ = beta_.^2 * sqrt(E*I/(rho*A));

% Use variable precision arithmetic for mode shapes to avoid numerical instability
w_hat_ = @(x,n) double( ...
    vpa(sin(beta_(n)*x)) - vpa(sinh(beta_(n)*x)) - ...
    ( (vpa(sin(beta_(n)*L)) - vpa(sinh(beta_(n)*L))) / ...
    (vpa(cos(beta_(n)*L)) - vpa(cosh(beta_(n)*L))) ) * ...
    (vpa(cos(beta_(n)*x)) - vpa(cosh(beta_(n)*x))) );

% Plot the vibration modes
x = linspace(0, L, 1000);

q4_modes = figure(1);clf;
t4 = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:4
    nexttile;
    hold on;
    for j = 1:3
        mode_idx = (i-1)*3 + j;
        if mode_idx <= N
            plot(x, w_hat_(x, mode_idx), 'DisplayName', sprintf('Mode %d', mode_idx), 'LineWidth', 1.5);
        end
    end
    hold off;
    grid on;
    xlabel('Position x [m]');
    ylabel('Displacement w(x) [m]');
    legend('Location', 'southwest');
end
if savefigs
    uniformFigureStyle(q4_modes, 'Q4_Beam_Modes', 20, 3/4);
end

% Display frequency and natural frequency results
fprintf('Mode\tBeta_n [1/m]\tOmega_n [rad/s]\n');
for n = 1:N
    fprintf('%d\t%.4f\t\t%.4f\n', n, beta_(n), omega_(n));
end

%% Question 5: Initial Deflection Simulation
% Define an arbitrary, but physically realistic initial deflection of the beam
% that complies to the double-clamped boundary conditions. Assume zero initial
% velocity and simulate for various values of N the deflection over a time
% period

% Initial deflection function
w0 = @(x) (1 * ((x/L).*(1 - x/L)).^3 .* sin(3*pi*x/L));

t_end = 0.2; % [s]
dt = 0.0002; % [s]
t = 0:dt:t_end;

x = linspace(0, L, 500); % [m]
N_values = [1, 3, 5];

w_xt = zeros(length(x), length(t),length(N_values));

for idx = 1:length(N_values)
    N = N_values(idx);

    % Calculate modal coefficients
    C_n = zeros(N, 1);
    for n = 1:N
        integrand = @(x) w0(x) .* w_hat_(x, n);
        C_n(n) = (1 / (rho * A * L)) * integral(integrand, 0, L);
    end
    % Simulate deflection over time

    for n = 1:N
        % Add the n-th modal contribution to the time-domain response.
        % For zero initial velocity, the modal expansion loses the sine term.
        w_xt(:, :, idx) = w_xt(:, :, idx) + C_n(n) * w_hat_(x, n)' * cos(omega_(n) * t);
    end
end

% Plot results over time
q5_deflection = figure(2);clf;
for idx = 1:length(N_values)
    N = N_values(idx);

    subplot(length(N_values), 1, idx);
    imagesc(t, x, w_xt(:, :, idx));
    axis xy;
    cb = colorbar;
    ylabel(cb, 'Deflection w(x,t) [m]','Rotation',270,'FontSize',10,'interpreter','latex');
    xlabel('Time [s]');
    ylabel('Position x [m]');
    title(sprintf('N=%d', N));
end

% Compare initial deflection to modal approximation
q5_initial = figure(3);clf;
hold on;
plot(x, w0(x), 'k--', 'DisplayName', 'Initial Deflection w_0(x)', 'LineWidth', 1.5);
for idx = 1:length(N_values)
    plot(x, w_xt(:,1,idx), 'DisplayName', sprintf('Modal Approx. N=%d', ...
        N_values(idx)), 'LineWidth', 1.5);
end
hold off;
grid on;
xlabel('Position x [m]');
ylabel('Deflection w(x,0) [m]');
legend('Location', 'southwest');

if savefigs
    uniformFigureStyle(q5_deflection, 'Q5_Deflection_Simulation', 18, 3/4);
    uniformFigureStyle(q5_initial, 'Q5_Initial_Comparison', 15, 1/2);
end

%% Question 8: Simulation Experiments with State-Space Model
% Use the state-space model from Question 7 to simulate beam deflections.
% Experiment with different inputs, approximation orders N, and initial conditions.
N_max_q8 = 20;

% Compute additional beta and omega values if needed
if length(beta_) < N_max_q8
    for i = length(beta_)+1:N_max_q8
        beta_(i) = double(vpasolve(eqn_beta, beta, (i+1)*pi/L));
    end
    omega_ = beta_.^2 * sqrt(E*I/(rho*A));

    % Redefine w_hat_ with expanded beta_ array
    w_hat_ = @(x,n) double( ...
        vpa(sin(beta_(n)*x)) - vpa(sinh(beta_(n)*x)) - ...
        ( (vpa(sin(beta_(n)*L)) - vpa(sinh(beta_(n)*L))) / ...
        (vpa(cos(beta_(n)*L)) - vpa(cosh(beta_(n)*L))) ) * ...
        (vpa(cos(beta_(n)*x)) - vpa(cosh(beta_(n)*x))) );
end

% Compute normalization constants for mode shapes
norm_const_q8 = zeros(N_max_q8, 1);
for n = 1:N_max_q8
    integrand = @(x) w_hat_(x, n).^2;
    norm_const_q8(n) = sqrt(integral(integrand, 0, L));
end

% Normalized mode shape function
phi_q8 = @(x,n) w_hat_(x,n) / norm_const_q8(n);

% Compute modal input coefficients b_n for state-space model
ell = @(x) double((x >= L_0 - W/2) & (x <= L_0 + W/2));
b_n = zeros(N_max_q8, 1);
for n = 1:N_max_q8
    integrand = @(x) ell(x) .* phi_q8(x, n);
    b_n(n) = (1/(rho*A)) * integral(integrand, 0, L);
end

% --- Experiment 1: Step Input Response ---
u_step = @(t) 100; % [N]
t_span = [0, 1]; % [s]
N_values_q8 = [2, 3, 4, 5, 10];
x_mid = L/2;

q8_step = figure(4); clf;
hold on;
colors_q8 = lines(length(N_values_q8));

for idx = 1:length(N_values_q8)
    N_q8 = N_values_q8(idx);
    fprintf('  N = %d modes...', N_q8);
    tic;

    % Build state-space matrices
    Omega_sq = diag(omega_(1:N_q8).^2);
    A_ss = [zeros(N_q8), eye(N_q8); -Omega_sq, zeros(N_q8)];
    B_ss = [zeros(N_q8, 1); b_n(1:N_q8)];

    % Initial state (zero deflection, zero velocity)
    x0_ss = zeros(2*N_q8, 1);

    % ODE function
    odefun = @(t, x) A_ss*x + B_ss*u_step(t);

    % Solve ODE
    [t_sol, x_sol] = ode45(odefun, t_span, x0_ss);
    a_coeff = x_sol(:, 1:N_q8);

    % Compute deflection at midpoint
    w_mid = zeros(length(t_sol), 1);
    for i = 1:length(t_sol)
        for n = 1:N_q8
            w_mid(i) = w_mid(i) + a_coeff(i, n) * phi_q8(x_mid, n);
        end
    end

    fprintf(' %.2f sec\n', toc);

    plot(t_sol, w_mid*1e3, 'LineWidth', 2, 'Color', colors_q8(idx,:), ...
        'DisplayName', sprintf('N=%d', N_q8));
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
tic;

N_q8 = 10;
omega_resonance = omega_(3);
u_sin = @(t) 50 * sin(omega_resonance * t); % [N]

% Build state-space matrices
Omega_sq = diag(omega_(1:N_q8).^2);
A_ss = [zeros(N_q8), eye(N_q8); -Omega_sq, zeros(N_q8)];
B_ss = [zeros(N_q8, 1); b_n(1:N_q8)];
x0_ss = zeros(2*N_q8, 1);

% Solve ODE
odefun = @(t, x) A_ss*x + B_ss*u_sin(t);
[t_sol, x_sol] = ode45(odefun, t_span, x0_ss);
a_coeff = x_sol(:, 1:N_q8);

% Spatial reconstruction at selected time instances
x_grid = linspace(0, L, 200);
t_snapshots = [0.25 0.50 0.75 1];

q8_resonance = figure(5); clf;
hold on;

for idx = 1:length(t_snapshots)
    [~, t_idx] = min(abs(t_sol - t_snapshots(idx)));

    % Reconstruct spatial deflection
    w_spatial = zeros(size(x_grid));
    for i = 1:length(x_grid)
        for n = 1:N_q8
            w_spatial(i) = w_spatial(i) + a_coeff(t_idx, n) * phi_q8(x_grid(i), n);
        end
    end

    plot(x_grid, w_spatial*1e3, 'LineWidth', 2, ...
        'DisplayName', sprintf('t=%.1fs', t_snapshots(idx)));
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
w0_q8 = @(x) sin(3*pi*x/L) .* (x/L).^3 .* (1 - x/L).^3;
u_zero = @(t) 0;
t_span_q8 = [0, 1];
N_values_q8_2 = [2, 3, 4, 5, 10];

q8_initial = figure(6); clf; hold on;

for idx = 1:length(N_values_q8_2)
    N_q8 = N_values_q8_2(idx);
    fprintf('  N = %d modes...', N_q8);
    tic;

    % Build state-space matrices
    Omega_sq = diag(omega_(1:N_q8).^2);
    A_ss = [zeros(N_q8), eye(N_q8); -Omega_sq, zeros(N_q8)];
    B_ss = [zeros(N_q8, 1); b_n(1:N_q8)];

    % Compute initial modal coefficients
    a0_q8 = zeros(N_q8, 1);
    for n = 1:N_q8
        integrand = @(x) w0_q8(x) .* phi_q8(x, n);
        a0_q8(n) = integral(integrand, 0, L);
    end

    x0_ss = [a0_q8; zeros(N_q8, 1)];

    % Solve ODE
    odefun = @(t, x) A_ss*x;
    [t_sol, x_sol] = ode45(odefun, t_span_q8, x0_ss);
    a_coeff = x_sol(:, 1:N_q8);

    % Compute deflection at midpoint
    w_mid = zeros(length(t_sol), 1);
    for i = 1:length(t_sol)
        for n = 1:N_q8
            w_mid(i) = w_mid(i) + a_coeff(i, n) * phi_q8(x_mid, n);
        end
    end

    fprintf(' %.2f sec\n', toc);

    plot(t_sol, w_mid*1e3, 'Color', colors_q8(idx,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('N=%d', N_q8));
    grid on;
    if idx == length(N_values_q8_2)
        grid on;
        ylabel('$w(L/2,t)$ [mm]', 'Interpreter', 'latex');
        xlabel('Time [s]', 'Interpreter', 'latex');
        legend('Location', 'southeast', 'Interpreter', 'latex');
    end
end

if savefigs
    uniformFigureStyle(q8_initial, 'Q8_Initial_Deflection', 18, 1/4);
end

%% Question 9: POD Basis Computation
% Generate snapshot data from a representative experiment with large N,
% compute POD basis via SVD, and implement POD-based state-space model.

% --- Step 1: Generate snapshot data with large N (reference simulation) ---
N_ref = N_max_q8;  % Use all available modes for reference
x_pod = linspace(0, L, 100)';  % Spatial grid (100 points)
t_pod = linspace(0, 1, 500);   % Time grid (500 snapshots)

% Build state-space matrices for reference model
Omega_sq_ref = diag(omega_(1:N_ref).^2);
A_ref = [zeros(N_ref), eye(N_ref); -Omega_sq_ref, zeros(N_ref)];
B_ref = [zeros(N_ref, 1); b_n(1:N_ref)];

% Combined input: step + sinusoidal at omega_1 to excite multiple modes
u_combined = @(t) 50 + 30*sin(omega_(1)*t) + 20*sin(omega_(3)*t);

% Initial deflection (same as Q8 Experiment 3)
a0_ref = zeros(N_ref, 1);
for n = 1:N_ref
    integrand = @(x) w0_q8(x) .* phi_q8(x, n);
    a0_ref(n) = integral(integrand, 0, L);
end
x0_ref = [a0_ref; zeros(N_ref, 1)];

% Solve ODE for reference solution
odefun_ref = @(t, x) A_ref*x + B_ref*u_combined(t);
[t_ref, x_ref] = ode45(odefun_ref, t_pod, x0_ref);
a_ref = x_ref(:, 1:N_ref);  % Modal coefficients [M_t x N_ref]

% Build snapshot matrix W: each column is w(x, t_k)
% W(i,k) = w(x_i, t_k) = sum_n a_n(t_k) * phi_n(x_i)
M_x = length(x_pod);
M_t = length(t_ref);

% Pre-compute mode shape matrix Phi_ref(i,n) = phi_n(x_i)
Phi_ref = zeros(M_x, N_ref);
for n = 1:N_ref
    Phi_ref(:, n) = phi_q8(x_pod, n);
end

% Snapshot matrix: W = Phi_ref * a_ref'
W = Phi_ref * a_ref';  % [M_x x M_t]

% --- Step 2: Compute POD basis via SVD ---
% SVD of snapshot matrix: W = U * S * V'
% POD modes are columns of U
[U_pod, S_pod, V_pod] = svd(W);

% Singular values (energy content)
sigma = diag(S_pod);
energy_content = cumsum(sigma.^2) / sum(sigma.^2);

% Display energy capture for different R values
fprintf('  Singular value energy capture:\n');
R_values_display = [1, 2, 3, 5, 10, 15, 20];
for R = R_values_display
    if R <= length(sigma)
        fprintf('    R=%2d: %.4f%% energy\n', R, energy_content(R)*100);
    end
end

% --- Plot: Singular value spectrum ---
q9_spectrum = figure(7); clf;
subplot(1,2,1);
h = stem(1:min(20, length(sigma)), sigma(1:min(20, length(sigma))), 'Color', colors_q8(1,:), 'LineWidth', 1.5);
% set(gca,'yscal','log')
grid on;
xlabel('Mode index $r$', 'Interpreter', 'latex');
ylabel('Singular value $\sigma_r$', 'Interpreter', 'latex');
title('POD Singular Value Spectrum', 'Interpreter', 'latex');

subplot(1,2,2);
plot(1:min(20, length(energy_content)), energy_content(1:min(20, length(energy_content)))*100, 'Color', colors_q8(1,:), 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
hold on;
yline(99, 'k--', 'LineWidth', 1);
yline(99.9, 'k:', 'LineWidth', 1);
hold off;
grid on;
xlabel('POD order $R$', 'Interpreter', 'latex');
ylabel('Cumulative energy [\%]', 'Interpreter', 'latex');
title('Energy Capture vs POD Order', 'Interpreter', 'latex');
legend({'Energy', '99\%', '99.9\%'}, 'Location', 'southeast', 'Interpreter', 'latex');

if savefigs
    uniformFigureStyle(q9_spectrum, 'Q9_POD_Spectrum', 18, 1/4);
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