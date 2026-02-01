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