
% Ryan Wood 4700 laser propagation model

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

InputParasL.E0 = 0e5;        % Input waveform parameters
InputParasL.we = 0;          % rotation frequency
InputParasL.t0 = 30e-13;      % time delay
InputParasL.wg = 10e-13;      % gaussian width
InputParasL.phi = 0;         % phase
InputParasL.rep = 500e-12;   % repeating pulse
InputParasR = 0;             % 

% Behavioural constants
n_g = 3.5;                   % index of refration of waveguide
vg = c_c/n_g*1e2;            % TWM cm/s group velocity
lambda = 1550e-9;            % wavelength
f0 = c_c/lambda;             % characteristic frequency
Ntr = 1e18;                  % for carrier frequency
g_fwhm = 3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 1e12;                  % resonant frequency

LGain = 0.05;                 % resonant gain
beta_spe = .3e-5;               % spontaneous emission beta
gamma = 1.0;                    % gamma SPE
SPE = 1;                        % spontaneous emission
taun = 1e-9;                    % spontaneous emission term tau_n
Zg = sqrt(c_mu_0/c_eps_0)/n_g;  % 
EtoP = 1/(Zg*f0*vg*1e-2*c_hb);  % conversion factor from electric field to photon density using Zg

% Waveguide and plotting parameters
plotN = 100;
L = 500e-4;             % cm length of waveguide
XL = [0,L];
YL = [-0*InputParasL.E0,300*InputParasL.E0];
RL = 0.5i;                   % Reflected wave left
RR = 0.5i;                   % Reflected wave right

Nz = 51;                    % number of z units in waveguide
dz = L/(Nz-1);               % cm unit of length
dt = dz/vg;                  % s unit of time
fsync = dt*vg/dz;            % Hz Synchronous frequency
Nt = floor(400*Nz);          % number of time units
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

% Grating
% kappa = grating(0,0,'');
kappa = 0*ones(size(z));

% Carrier Equation
N = ones(size(z))*Ntr;          %
Nave = nan(1,Nt);
Nave(1) = mean(N);              % average of N
gain = vg*2.5e-16;              % G0 gain
eVol = 1.5e-10*c_q;             % volume of electrons
Ion = (dt*Nt)*1/5;                  % current is added between Ion and Ioff
Ioff = (dt*Nt);                    % 
I_off = 0.024;                  % the value of current provided when turned off
I_on = 0.15;                     % vs. when turned on

% Initial Conditions
Ef1 = @SourceFct;            % forward envelope child of sourcefct
ErN = @SourceFct;            % reflected envelope child of sourcefct

t = 0;                       % set time to 0
time(1) = t;                 % set first index of time array to t

InputL(1) = Ef1(t,InputParasL);     % first index of left input array set to gaussian from sourcefct
InputR(1) = ErN(t,InputParasR);     % first index of right input array set to gaussian from sourcefct

OutputR(1) = Ef(Nz);                % first index of right output set to value of reflected electric field envelope on right end of waveguide
OutputL(1) = Er(1);                 % first index of left output set to the value of the incident electric field envelope on left end of waveguide

Ef(1) = InputL(1);                  % set forward envelope to left input initial condition
Er(Nz) = InputR(1);                 % set reverse envelope to right input initial condition

% Plotting
figure('name','Fields')             % set up plotting
subplot(3,2,1)
plot(z*10000,real(Ef),'r');hold on         % plot real part of forward envelope
xlabel('z(\mum)')
ylabel('E_f')
hold off

subplot(3,2,2)
plot(z*10000,N,'r');     hold on           % plot N
xlabel('z(\mum)')
ylabel('N')
hold off

subplot(3,2,[3,4])                  % plot Nave
plot(time*1e12,Nave,'b');  hold on    
xlabel('z(\mum)')
ylabel('Nave')
hold off

subplot(3,2,[5,6])                  % plot inputs and outputs
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR),'g');
plot(time*1e12,real(InputR),'b');
plot(time*1e12,real(OutputL),'p--');
xlabel('time(ps)')
ylabel('E')
hold off

% Main loop

for i = 2:Nt                    % loop while time is between 2 and Nt (number of time units)
    t = dt*(i-1);               % set time to next time unit
    time(i) = t;                % set each time array index to current time each loop

    A = sqrt(gamma*beta_spe*c_hb*f0*L*1e-2/taun)/(2*Nz);
    if SPE > 0
        Ff = (randn(Nz,1)+1i*randn(Nz,1))*A;
        Fr = (randn(Nz,1)+1i*randn(Nz,1))*A;
    else
        Ff = (ones(Nz,1))*A;
        Fr = (ones(Nz,1))*A;
    end

  

    % Setting up equations
    InputL(i) = Ef1(t,InputParasL);     % propagate envelope from left end of waveguide
    InputR(i) = ErN(t,0);

    Ef(1) = InputL(i) + RL*Er(1);       % first index of forward envelope is current left input plus the left reflected envelope
    Er(Nz) = InputR(i) + RR*Ef(Nz);     % last index of reverse envelope is current right input plus the right reflected envelope

    Pf(1) = 0;  
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;

    S = (abs(Ef).^2 +abs(Er).^2).*EtoP*1e-6;            % value of S for stimulated emission
    if t < Ion || t > Ioff                              % decide when to turn the current on or off based on time
        I_injv = I_off;                                 % turn off
    else                                                %
        I_injv = I_on;                                  % turn on
    end                                                 %      

    Stim = gain.*(N-Ntr).*S;                            % stimulated emission term
    N = (N + dt*(I_injv/eVol - Stim))./(1+dt/taun);     % N(z) update equation
    Nave(i) = mean(N);                                  % record average value of N each step

  % Exponential spatial growth term and stimulated emission
    alpha = 0;                                  % loss term
    beta_r = 0;                                 % real wave beta
    gain_z = gain.*(N - Ntr)./vg;                % gain term due to stimulated emission
    beta_i = (gain_z-alpha)./2;                 % imaginary wave beta with gain/loss terms
    beta = 1i*beta_i;                           % set length of waveguide to beta
    exp_det = exp(-1i*dz*beta);                 % exponential growth term along z, rotating with length

    EsF = Ff*abs(SPE).*sqrt(N.*1e6);
    EsR = Fr*abs(SPE).*sqrt(N.*1e6);

    % Update equations
    Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);             % forward polarization term using trapezoidal rule
    Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);          
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);             % reverse polarization term calculated with trapezoidal rule


    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain*(Ef(2:Nz-1)-Pf(2:Nz-1));            % forward envelope with polarization and gain terms
    Er(2:Nz-1) = Er(2:Nz-1) - LGain*(Er(2:Nz-1)-Pr(2:Nz-1));            % reverse envelope with polarization and gain terms

    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Efp(1:Nz-1) + 1i*dz*kappa(1:Nz-1).*Er(1:Nz-1);        % updated value of forward envelope
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(2:Nz).*Efp(2:Nz);        % updated value of reverse envelope
    
    OutputR(i) = Ef(Nz)*(1-RR);         % set right output to the last index of forward envelope scaled by 1 minus the reflection coefficient
    OutputL(i) = Er(1)*(1-RL);          % set left output to the first index of reverse envelope scaled by 1 minutes the reflection coefficient

    Ef = Ef + EsF;
    Er = Er + EsR;

    % Reset to previous values to avoid instability
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;


    % Live plotting of input and output signals
    if mod(i,plotN) == 0

        gain_z

        % Forward real and imaginary waveform
        subplot(3,2,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--');
        xlim(XL*1e4)
        ylim([-4e6,4e6])
        xlabel('z(\mum)')
        ylabel('E_f')
        % legend('\Re','\Im')
        hold off

        % Value of N
        subplot(3,2,2)
        plot(z*10000,N,'r'); hold on
        xlim(XL*1e4)
        ylim([0,5e18])
        xlabel('z(\mum)')
        ylabel('N')
        % legend('\Re','\Im')
        hold off

        % Reverse real and imaginary waveform
        subplot(3,2,[3,4])
        plot(time*1e12,Nave,'b'); hold on
        xlim([0,Nt*dt*1e12])
        ylim([0,5e18])
        xlabel('time(ps)')
        ylabel('Nave')
        % legend('\Re','\Im')
        hold off

        % Live plotting of left and right-side I/O waveforms
        subplot(3,2,[5,6])
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12,real(OutputR),'g');
        plot(time*1e12,real(InputR),'b');
        plot(time*1e12,real(OutputL),'m');
        xlim([0,Nt*dt*1e12])
        ylim([-5e6, 5e6])
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input','Right Output','Right Input','Left Output')
        hold off

        pause(0.01)
    end
end

% Frequency domain analysis of signal

omega = fftshift(wspace(time));   % Find phase shift of the fourier transform

% OutputR = OutputR - mean(OutputR);
% len = length(OutputR);
% Fs = 1/0.0117;                      %sampling
% Fn = Fs/2;
% fftOut = fftshift(fft(OutputR));
% Iv = 1:length(omega);
% fftOut_a = double(abs(fftOut(Iv))*2);
% fftOutput = sgolayfilt(fftOut_a,2,21); %savitsky-golay filtering method of FFT




fftOutput_raw = fftshift(fft(OutputR));
fftOutput = envelope(abs(fftOutput_raw), 2, 'peak');  %envelope filtering method of FFT
[yupperFFT,ylowerFFT] = envelope(abs(fftOutput_raw), 2, 'peak');


[yupper,ylower] = envelope(real(OutputR), 300, 'peak');             %startup transient of laser - ylower



figure('name','Frequency Analysis')

% Plot the input and output waveforms over time
subplot(1,3,1)
plot(time*1e12,real(InputL),'g'); hold on
plot(time*1e12,ylower,'r');
xlim([0,Nt*dt*1e12])
ylim([min(ylower), max(ylower)])
xlabel('time(ps)')
ylabel('Right Output')
legend('Input','Output')
hold off

% Plot the fourier transform of the input and output waveforms
% nonZero_in = fftInput ~= 0;
% db_in(nonZero_in) = 20*log10(abs(fftInput(nonZero_in)));
nonZero_out = fftOutput ~= 0;
db_out(nonZero_out) = 20*log10(abs(fftOutput(nonZero_out)));

subplot(1,3,2)
% plot(omega,db_in,'b'); hold on
plot(omega*1e-12,20*log10(fftOutput),'g'); hold off
xlim([-5,5])
% ylim([min(fftOutput),max(fftOutput)])
xlabel('GHz')
ylabel('|E| (dB)')
legend('Input','Output')

% Plot the phase shift of both waveforms
fftOutput = fftshift(fft(abs(OutputR)));     % FFT of output
fftInput = fftshift(fft(InputL));       % FFT of input

% yMin = min(min(unwrap(angle(fftOutput))), min(unwrap(angle(fftInput))));
% yMax = max(max(unwrap(angle(fftOutput))),max(unwrap(angle(fftInput))));

subplot(1,3,3)
% plot(omega,unwrap(angle(fftInput)),'b'); hold on
plot(omega*1e-12,unwrap(angle(fftOutput)),'b'); hold off
xlim([-5,5])
% ylim([yMin,yMax])
xlabel('GHz')
ylabel('phase (E)')
legend('Input','Output')


























