% Ryan Wood 4700 laser propagation model

set(0,'defaultaxesfontsize',20)
set(0,'DefaultFigureWindowStyle','docked')
set(0,"DefaultLineLineWidth",2)
set(0,'Defaultaxeslinewidth',2)

set(0,'DefaultFigureWindowStyle','docked')

c_c = 299792458;             % m/s TWM speed of light
c_eps_0 = 8.8542149e-12;     % F/m vacuum permittivity
c_eps_0_cm = c_eps_0/100;    % F/cm
c_mu_0 = 1/c_eps_0/c_c^2;    % H/m
c_q = 1.60217653e-19;        % Coulombs
c_hb = 1.05457266913e-34;    % h bar
c_h = c_hb*2*pi;             % h

InputParasL.E0 = 1e5;        % Input waveform parameters
InputParasL.we = 0;          % imaginary rotation frequency
InputParasL.t0 = 2e-12;      % time delay
InputParasL.wg = 5e-13;      % gaussian width
InputParasL.phi = 0;         % phase
InputParasR = 0;             % 


n_g = 3.5;                   % index of refration of waveguide
vg = c_c/n_g*1e2;            % TWM cm/s group velocity
lambda = 1550e-9;

plotN = 10;

L = 1000e-6*1e2;             % cm length of waveguide
XL = [0,L];
YL = [0,InputParasL.E0];

RL = 0.9i;                   % Reflected wave left
RR = 0.9i;                   % Reflected wave right

Nz = 500;                    % number of z units in waveguide
dz = L/(Nz-1);               % cm unit of length
dt = dz/vg;                  % s unit of time
fsync = dt*vg/dz;            % Hz Synchronous frequency

Nt = floor(2*Nz);            % number of units of time is twice that of length
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

figure('name','Fields')             % set up plotting
subplot(3,1,1)
plot(z*10000,real(Ef),'r');         % plot real part of forward envelope
hold off
xlabel('z(\mum)')
ylabel('E_f')
subplot(3,1,2)
plot(z*10000,real(Er),'b');         % plot real part of reverse envelope
xlabel('z(\mum)')
ylabel('E_r')
hold off
subplot(3,1,3)
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR),'r--');
plot(time*1e12,real(InputR),'b'); hold on
plot(time*1e12,real(OutputL),'b--');
xlabel('time(ps)')
ylabel('E')

hold off

for i = 2:Nt                    % loop while time is between 2 and Nt (number of time units)
    t = dt*(i-1);               % set time to next time unit
    time(i) = t;                % set each time array index to current time each loop

    InputL(i) = Ef1(t,InputParasL);     % propagate envelope from left end of waveguide
    InputR(i) = ErN(t,0);

    Ef(1) = InputL(i) + RL*Er(1);       % first index of forward envelope is current left input plus the left reflected envelope
    Er(Nz) = InputR(i) + RR*Ef(Nz);     % last index of reverse envelope is current right input plus the right reflected envelope

    Ef(2:Nz) = fsync*Ef(1:Nz-1);        % value of forward envelope past the first index
    Er(1:Nz-1) = fsync*Er(2:Nz);        % value of the reverse envelope previous to the final index

    OutputR(i) = Ef(Nz)*(1-RR);         % set right output to the last index of forward envelope scaled by 1 minus the reflection coefficient
    OutputL(i) = Er(1)*(1-RL);          % set left output to the first index of reverse envelope scaled by 1 minutes the reflection coefficient

    if mod(i,plotN) == 0
        subplot(3,1,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re','\Im')
        hold off
        subplot(3,1,2)
        plot(z*10000,real(Er),'b'); hold on
        plot(z*10000,imag(Er),'b--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re','\Im')

        hold off
        subplot(3,1,3)
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12,real(OutputR),'g');
        plot(time*1e12,real(InputR),'b');
        plot(time*1e12,real(OutputL),'m');
        xlim([0,Nt*dt*1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input','Right Output','Right Input','Left Output', 'Location','East')
        hold off
        pause(0.01)
    end
end












