% Ryan Wood 4700 laser propagation model - Modified to 1D DFB Laser

clear all;
close all;
set(0,'defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','docked')
set(0,"DefaultLineLineWidth",2)
set(0,'Defaultaxeslinewidth',2)
set(0,'DefaultFigureWindowStyle','docked')

% Constants and terms
c_c = 299792458;             % m/s speed of light
c_eps_0 = 8.8542149e-12;     % F/m vacuum permittivity
c_eps_0_cm = c_eps_0/100;    % F/cm
c_mu_0 = 1/c_eps_0/c_c^2;    % H/m
c_q = 1.60217653e-19;        % Coulombs
c_hb = 1.05457266913e-34;    % h bar
c_h = c_hb*2*pi;             % h

InputParasL.E0 = 0e5;        % Input waveform parameters
InputParasL.we = 0;          % rotation frequency
InputParasL.t0 = 30e-13;     % time delay
InputParasL.wg = 10e-13;     % gaussian width
InputParasL.phi = 0;         % phase
InputParasL.rep = 500e-12;   % repeating pulse
InputParasR = 0;             % 

% Behavioural constants
n_g = 3.5;                   % index of refraction of waveguide
vg = c_c/n_g*1e2;            % cm/s group velocity
lambda = 1550e-9;            % wavelength
f0 = c_c/lambda;             % characteristic frequency
Ntr = 1e18;                  % carrier density transparency
g_fwhm = 3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 1e12;                  % resonant frequency

LGain = 0.05;                % resonant gain
beta_spe = .3e-5;            % spontaneous emission beta
gamma = 1.0;                 % gamma SPE
SPE = 10;                    % spontaneous emission
taun = 1e-9;                 % spontaneous emission term tau_n
Zg = sqrt(c_mu_0/c_eps_0)/n_g;  
EtoP = 1/(Zg*f0*vg*1e-2*c_hb);  % conversion factor from E-field to photon density

% Waveguide and plotting parameters
plotN = 100;
L = 100e-4;                  % cm length of waveguide
XL = [0,L];
YL = [-0*InputParasL.E0,300*InputParasL.E0];
RL = 0;%.9i;                      % Reflected wave left (set to 0 for DFB)
RR = 0;%.9i;                      % Reflected wave right (set to 0 for DFB)

Nz = 51;                     % number of z units in waveguide
dz = L/(Nz-1);               % cm unit of length
dt = dz/vg;                  % s unit of time
fsync = dt*vg/dz;            % Hz Synchronous frequency
Nt = floor(600*Nz);          % number of time units
tmax = Nt*dt;                % time length of the simulation
t_L = dt*Nz;                 % time to travel length

z = linspace(0,L,Nz).';      % Nz points, Nz-1 segments
time = nan(1,Nt);            % 1xNt array of NaN
InputL = nan(1,Nt);          %
InputR = nan(1,Nt);          %
OutputL = nan(1,Nt);         %
OutputR = nan(1,Nt);         %

Ef = zeros(size(z));         % envelope of forward electric field
Er = zeros(size(z));         % envelope of reverse electric field
Pf = zeros(size(z));         % forward polarization wave
Pr = zeros(size(z));         % reverse polarization wave

Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;

% DFB Grating
kappa_0 = 1300;               % Coupling coefficient (cm^-1), adjust as needed
Lambda = lambda/(2*n_g);         % Grating period (cm), tied to Bragg wavelength
kappa = kappa_0 * cos(2 * pi * z / Lambda);  % Sinusoidal grating

% Carrier Equation
N = ones(size(z))*Ntr;       %
Nave = nan(1,Nt);
Nave(1) = mean(N);           % average of N
gain = vg*2.5e-16;           % G0 gain
eVol = 1.5e-10*c_q;          % volume of electrons
Ion = (dt*Nt)*1/7;           % current is added between Ion and Ioff
Ioff = (dt*Nt);              % 
I_off = 0.024;               % current off value
I_on = 0.15;                  % current on value

% Initial Conditions
Ef1 = @SourceFct;            % forward envelope child of sourcefct
ErN = @SourceFct;            % reflected envelope child of sourcefct

t = 0;                       % set time to 0
time(1) = t;                 % set first index of time array to t

InputL(1) = Ef1(t,InputParasL);     % initial left input
InputR(1) = ErN(t,InputParasR);     % initial right input
OutputR(1) = Ef(Nz);                % initial right output
OutputL(1) = Er(1);                 % initial left output

Ef(1) = InputL(1);                  % set forward envelope to left input
Er(Nz) = InputR(1);                 % set reverse envelope to right input

% Plotting setup (unchanged)
figure('name','Fields')
subplot(3,2,1)
plot(z*10000,real(Ef),'r'); hold on
xlabel('z(\mum)')
ylabel('E_f')
hold off

subplot(3,2,2)
plot(z*10000,N,'r'); hold on
xlabel('z(\mum)')
ylabel('N')
hold off

subplot(3,2,[3,4])
plot(time*1e12,Nave,'b'); hold on    
xlabel('z(\mum)')
ylabel('N_a_v_e')
hold off

subplot(3,2,[5,6])
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR),'g');
plot(time*1e12,real(InputR),'b');
plot(time*1e12,real(OutputL),'p--');
xlabel('time(ps)')
ylabel('E')
hold off

% Main loop
for i = 2:Nt
    t = dt*(i-1);
    time(i) = t;

    A = sqrt(gamma*beta_spe*c_hb*f0*L*1e-2/taun)/(2*Nz);
    if SPE > 0
        Ff = (randn(Nz,1)+1i*randn(Nz,1))*A;
        Fr = (randn(Nz,1)+1i*randn(Nz,1))*A;
    else
        Ff = (ones(Nz,1))*A;
        Fr = (ones(Nz,1))*A;
    end

    % Setting up equations
    InputL(i) = Ef1(t,InputParasL);
    InputR(i) = ErN(t,0);

    Ef(1) = InputL(i) + RL*Er(1);       % Boundary condition
    Er(Nz) = InputR(i) + RR*Ef(Nz);     % Boundary condition

    Pf(1) = 0;  
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;

    S = (abs(Ef).^2 + abs(Er).^2).*EtoP*1e-6;  % Stimulated emission source
    if t < Ion || t > Ioff
        I_injv = I_off;
    else
        I_injv = I_on;
    end

    Stim = gain.*(N-Ntr).*S;
    N = (N + dt*(I_injv/eVol - Stim))./(1+dt/taun);
    Nave(i) = mean(N);

    % Exponential spatial growth term and stimulated emission
    alpha = 0;                          % loss term
    gain_z = gain.*(N - Ntr)./vg;       % gain term
    beta_i = (gain_z-alpha)./2;         % imaginary wave beta
    beta = 1i*beta_i;
    exp_det = exp(-1i*dz*beta);

    EsF = Ff*abs(SPE).*sqrt(N.*1e6);
    EsR = Fr*abs(SPE).*sqrt(N.*1e6);

    % Update polarization terms
    Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
    Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);

    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain*(Ef(2:Nz-1)-Pf(2:Nz-1));
    Er(2:Nz-1) = Er(2:Nz-1) - LGain*(Er(2:Nz-1)-Pr(2:Nz-1));

    % DFB-specific update with grating coupling
    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Efp(1:Nz-1) + 1i*dz*kappa(1:Nz-1).*Er(1:Nz-1);
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(2:Nz).*Efp(2:Nz);

    OutputR(i) = Ef(Nz)*(1-RR);
    OutputL(i) = Er(1)*(1-RL);

    Ef = Ef + EsF;
    Er = Er + EsR;

    % Reset to previous values
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;

    % Live plotting (unchanged)
    if mod(i,plotN) == 0

        gain_z

        subplot(3,2,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--');
        xlim(XL*1e4)
        ylim([-4e6,4e6])
        xlabel('z(\mum)')
        ylabel('E_f')
        hold off

        subplot(3,2,2)
        plot(z*10000,N,'r'); hold on
        xlim(XL*1e4)
        % ylim([0,1.5e18])
        xlabel('z(\mum)')
        ylabel('N')
        hold off

        subplot(3,2,[3,4])
        plot(time*1e12,Nave,'b'); hold on
        xlim([0,Nt*dt*1e12])
        % ylim([0,1.5e18])
        xlabel('time(ps)')
        ylabel('N_a_v_e')
        hold off

        subplot(3,2,[5,6])
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12,real(OutputR),'g');
        plot(time*1e12,real(InputR),'b');
        plot(time*1e12,real(OutputL),'m');
        xlim([0,Nt*dt*1e12])
        ylim([-5e6, 5e6])
        xlabel('time(ps)')
        ylabel('E')
        hold off

        pause(0.01)
    end
end

% Frequency domain analysis
omega = fftshift(wspace(time));  % Frequency array (in Hz)
fftOutput_raw = fftshift(fft(OutputR));  % FFT of the right output
fftOutput = envelope(abs(fftOutput_raw), 30, 'peak');  % Envelope of the FFT
[yupperFFT, ylowerFFT] = envelope(abs(fftOutput_raw), 30, 'peak');

% Find the peak frequency
[~, peak_idx] = max(abs(fftOutput));  % Index of the peak amplitude
peak_freq = omega(peak_idx);  % Peak frequency in Hz

% Shift the frequency axis to center the peak at 0
omega_shifted = omega - peak_freq;  % Shifted frequency array (Hz)
omega_shifted_GHz = omega_shifted * 1e-12;  % Convert to GHz

% Compute dB scale for the FFT output
nonZero_out = fftOutput ~= 0;
db_out = zeros(size(fftOutput));
db_out(nonZero_out) = 20 * log10(abs(fftOutput(nonZero_out)));

% Envelope of the time-domain signal (unchanged)
[yupper, ylower] = envelope(real(OutputR), 300, 'peak');

% Plotting
figure('name', 'Frequency Analysis')

% Plot the input and output waveforms over time (unchanged)
subplot(1, 3, 1)
plot(time*1e12, real(InputL), 'g'); hold on
plot(time*1e12, ylower, 'r');
xlim([0, Nt*dt*1e12])
ylim([min(ylower), max(ylower)])
xlabel('time (ps)')
ylabel('Right Output')
legend('Input', 'Output')
hold off

% Plot the Fourier transform centered at the peak frequency
subplot(1, 3, 2)
plot(omega_shifted_GHz, db_out, 'g'); hold off
xlim([-f0*2e-15, f0*2e-15])  % Display ±5 GHz around the peak
xlabel('Frequency Offset (GHz)')
ylabel('20log_1_0(|E|)')
legend('Output')

% Plot the phase of the FFT centered at the peak frequency
subplot(1, 3, 3)
fftOutput = fftshift(fft(abs(OutputR)));  % FFT of the absolute output
plot(omega_shifted_GHz, unwrap(angle(fftOutput)), 'b'); hold off
xlim([-5, 5])  % Display ±5 GHz around the peak
xlabel('Frequency Offset (GHz)')
ylabel('Phase (E)')
legend('Output')