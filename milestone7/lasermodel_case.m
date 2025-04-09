% Ryan Wood 4700 laser propagation model - Modified for 4 specific cases

clear all;
close all;
set(0,'defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','docked')
set(0,"DefaultLineLineWidth",2)
set(0,'Defaultaxeslinewidth',2)
set(0,'DefaultFigureWindowStyle','docked')

% Constants and terms
c_c = 299792458;             % m/s TWM speed of light
c_eps_0 = 8.8542149e-12;     % F/m vacuum permittivity
c_eps_0_cm = c_eps_0/100;    % F/cm
c_mu_0 = 1/c_eps_0/c_c^2;    % H/m
c_q = 1.60217653e-19;        % Coulombs
c_hb = 1.05457266913e-34;    % h bar
c_h = c_hb*2*pi;             % h

% Input waveform parameters (base case)
InputParasL.E0 = 1e5;        % Amplitude (will be overwritten)
InputParasL.we = 0;          % rotation frequency
InputParasL.t0 = 30e-13;     % time delay
InputParasL.wg = 10e-13;     % gaussian width
InputParasL.phi = 0;         % phase
InputParasL.rep = 500e-12;   % repeating pulse
InputParasR = 0;

% Behavioural constants
n_g = 3.5;                   % index of refraction of waveguide
vg = c_c/n_g*1e2;            % TWM cm/s group velocity
lambda = 1550e-9;            % wavelength
f0 = c_c/lambda;             % characteristic frequency
Ntr = 1e18;                  % for carrier frequency
g_fwhm = 3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 1e12;                  % resonant frequency
LGain = 0.05;                % resonant gain (increased to enable lasing)
beta_spe = .3e-5;            % spontaneous emission beta
gamma = 1.0;                 % gamma SPE
SPE = 10;                    % spontaneous emission (enabled)
taun = 1e-9;                 % spontaneous emission term tau_n
Zg = sqrt(c_mu_0/c_eps_0)/n_g;
EtoP = 1/(Zg*f0*vg*1e-2*c_hb);  % conversion factor from electric field to photon density

% Waveguide and plotting parameters
plotN = 100;
L = 1000e-6*1e2;             % cm length of waveguide
XL = [0,L];
YL = [-0*InputParasL.E0,300*InputParasL.E0];
RL = 0;                      % Reflected wave left
RR = 0;                      % Reflected wave right

Nz = 100;                    % number of z units in waveguide
dz = L/(Nz-1);               % cm unit of length
dt = dz/vg;                  % s unit of time
fsync = dt*vg/dz;            % Hz Synchronous frequency
Nt = floor(400*Nz);          % number of time units
tmax = Nt*dt;                % time length of the simulation
t_L = dt*Nz;                 % time to travel length

z = linspace(0,L,Nz).';      % Nz points, Nz-1 segments

% Grating
kappa = 0*ones(size(z));

% Define the 4 specific cases
cases = [
    0.05, 0;      % Case 1: Low I_on, no input
    0.05, 2e5;    % Case 2: Low I_on, high input
    0.05, 5e5;      % Case 3: High I_on, no input
    0.10, 0;% Case 4: High I_on, high input
    0.10, 2e5;
    0.10, 5e5
    ];
num_cases = size(cases, 1);

% Storage for results
OutputR_all = cell(num_cases, 1);
OutputL_all = cell(num_cases, 1);
Nave_all = cell(num_cases, 1);
time_all = cell(num_cases, 1);
case_labels = cell(num_cases, 1);

% Run simulations for each case
for case_idx = 1:num_cases
    % Set parameters for this case
    I_on = cases(case_idx, 1);
    I_off = 0.024;  % Fixed I_off
    InputParasL.E0 = cases(case_idx, 2);

    % Reset arrays
    time = nan(1,Nt);
    InputL = nan(1,Nt);
    InputR = nan(1,Nt);
    OutputL = nan(1,Nt);
    OutputR = nan(1,Nt);
    Ef = zeros(size(z));
    Er = zeros(size(z));
    Pf = zeros(size(z));
    Pr = zeros(size(z));
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;
    N = ones(size(z))*Ntr;
    Nave = nan(1,Nt);
    Nave(1) = mean(N);
    gain = vg*2.5e-16;
    eVol = 1.5e-10*c_q;

    % Initial Conditions
    t = 0;
    time(1) = t;
    InputL(1) = SourceFct(t,InputParasL);
    InputR(1) = SourceFct(t,InputParasR);
    OutputR(1) = Ef(Nz);
    OutputL(1) = Er(1);
    Ef(1) = InputL(1);
    Er(Nz) = InputR(1);

    % Main loop
    for k = 2:Nt
        t = dt*(k-1);
        time(k) = t;

        A = sqrt(gamma*beta_spe*c_hb*f0*L*1e-2/taun)/(2*Nz);
        if SPE > 0
            Ff = (randn(Nz,1)+1i*randn(Nz,1))*A;
            Fr = (randn(Nz,1)+1i*randn(Nz,1))*A;
        else
            Ff = (ones(Nz,1))*A;
            Fr = (ones(Nz,1))*A;
        end

        InputL(k) = SourceFct(t,InputParasL);
        InputR(k) = SourceFct(t,InputParasR);

        Ef(1) = InputL(k) + RL*Er(1);
        Er(Nz) = InputR(k) + RR*Ef(Nz);

        Pf(1) = 0;
        Pf(Nz) = 0;
        Pr(1) = 0;
        Pr(Nz) = 0;
        Cw0 = -LGamma + 1i*Lw0;

        S = (abs(Ef).^2 + abs(Er).^2).*EtoP*1e-6;

        % Pulse stream for I_D
        pulse_period = 500e-12;
        half_period = pulse_period / 2;
        if mod(t, pulse_period) < half_period
            I_injv = I_on;
        else
            I_injv = I_off;
        end

        Stim = gain.*(N-Ntr).*S;
        N = (N + dt*(I_injv/eVol - Stim))./(1+dt/taun);
        Nave(k) = mean(N);

        alpha = 0;
        beta_r = 0;
        gain_z = gain.*(N - Ntr)./vg;
        beta_i = (gain_z-alpha)./2;
        beta = 1i*beta_i;
        exp_det = exp(-1i*dz*beta);

        EsF = Ff*abs(SPE).*sqrt(N.*1e6);
        EsR = Fr*abs(SPE).*sqrt(N.*1e6);

        Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
        Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
        Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
        Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);

        Ef(2:Nz-1) = Ef(2:Nz-1) - LGain*(Ef(2:Nz-1)-Pf(2:Nz-1));
        Er(2:Nz-1) = Er(2:Nz-1) - LGain*(Er(2:Nz-1)-Pr(2:Nz-1));

        Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Efp(1:Nz-1) + 1i*dz*kappa(1:Nz-1).*Er(1:Nz-1);
        Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(2:Nz).*Efp(2:Nz);

        OutputR(k) = Ef(Nz)*(1-RR);
        OutputL(k) = Er(1)*(1-RL);

        Ef = Ef + EsF;
        Er = Er + EsR;

        Efp = Ef;
        Erp = Er;
        Pfp = Pf;
        Prp = Pr;
    end

    % Store results
    OutputR_all{case_idx} = OutputR;
    OutputL_all{case_idx} = OutputL;
    Nave_all{case_idx} = Nave;
    time_all{case_idx} = time;
    case_labels{case_idx} = sprintf('Ion=%.2f, E0=%.1e', I_on, InputParasL.E0);
end

% Plotting Results
% figure('name', 'Effect of I_D and Input Pulse (4 Cases)')
% % Time Domain: OutputR
% subplot(2,2,1)
% for i = 1:num_cases
%     plot(time_all{i}*1e12, real(OutputR_all{i}), 'DisplayName', case_labels{i}); hold on
% end
% xlim([0, Nt*dt*1e12])
% xlabel('time (ps)')
% ylabel('OutputR (Real)')
% title('Output Power vs. I_D and E0')
% legend('show')
% hold off

% Time Domain: Nave
subplot(2,2,[1,2])
for i = 1:num_cases
    plot(time_all{i}*1e12, Nave_all{i}, 'DisplayName', case_labels{i}); hold on
end
xlim([0, Nt*dt*1e12])
xlabel('time (ps)')
ylabel('Nave')
title('Average Carrier Density vs. I_D and E0')
legend('show')
hold off

% Frequency Domain
subplot(2,2,[3,4])
for i = 1:num_cases
    omega = fftshift(wspace(time_all{i}));
    fftOutput = fftshift(fft(OutputR_all{i}));
    plot(omega*1e-12, 20*log10(abs(fftOutput)), 'DisplayName', case_labels{i}); hold on
end
xlim([-5, 5])  % ±5 GHz
xlabel('Frequency (GHz)')
ylabel('|E| (dB)')
title('Spectrum vs. I_D and E0')
legend('show')
hold off