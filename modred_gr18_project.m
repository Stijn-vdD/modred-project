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

savefigs = false;

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
    ylabel('Displacement w(x)');
    legend('Location', 'best');
end
if savefigs
    uniformFigureStyle(q4_modes, 'Beam_Modes');
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
w0 = @(x) (0.1 * ((x/L).*(1 - x/L)).^3 .* sin(3*pi*x/L)); % 0.1 m max deflection

t_end = 0.2; % [s]
dt = 0.0002; % [s]
t = 0:dt:t_end;

x = linspace(0, L, 500); % [m]
N_values = [1, 3, 5];

q5_deflection = figure(2);clf;
for idx = 1:length(N_values)
    N = N_values(idx);

    % Calculate modal coefficients
    C_n = zeros(N, 1);
    for n = 1:N
        integrand = @(x) w0(x) .* w_hat_(x, n);
        C_n(n) = (1 / (rho * A * L)) * integral(integrand, 0, L);
    end
    % Simulate deflection over time
    w_xt = zeros(length(x), length(t));
    w_xt(:, 1) = w0(x)'; % Initial deflection at t=0
    for n = 1:N
        % Add the n-th modal contribution to the time-domain response.
        % For zero initial velocity, the modal expansion loses the sine term.
        w_xt(:,2:end) = w_xt(:,2:end) + C_n(n) * w_hat_(x, n)' * cos(omega_(n) * t(2:end));
    end
    % Plot results
    subplot(length(N_values), 1, idx);
    imagesc(t, x, w_xt);
    axis xy;
    colorbar;
    xlabel('Time [s]');
    ylabel('Position x [m]');
    title(sprintf('N=%d', N));
end
if savefigs
    uniformFigureStyle(q5_deflection, 'Beam_Deflection_Simulation');
end

%% Figure export function
function result = uniformFigureStyle(figureHandle,fileName)
% Author: S. van den Dungen, 2025

% Change the figure size
figureHandle.Units               = 'centimeters';
figureHandle.Position(3)         = 20;
figureHandle.Position(4)         = (3/4)*figureHandle.Position(3);

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

    % Change font
    %     ax.FontName = 'Myriad Pro';
    ax.FontSize = 10;

    result = result + 1;
end

% 'Remove' all borders from legends
leg = findall(figureHandle, 'Tag', 'legend');
if numel(leg)
    % Set the edge color to white ,without removing its white
    % background
    set(leg, 'EdgeColor', [1 1 1]);
end
legend('Location','best')

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