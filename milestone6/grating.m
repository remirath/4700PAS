% Ryan Wood 101195542 milestone5
function g = grating(start,stop,type)

L = 1000e-6*1e2;             % cm length of waveguide
Nz = 500;                    % number of points along z axis
z = linspace(0,L,Nz).';      % Nz points, Nz-1 segments
kappa = zeros(size(z));      % kappa as matrix

% Uniform
if (strcmp(type,'uniform'))
    kappa0 = 100;                % coupling coefficient
    temp = linspace(0,5,(Nz*stop)-(Nz*start)).'; % range from 0 to 1 over grating length
    kappa(round(Nz*start):round(Nz*start)+length(temp)-1) = 1;
    kappa = kappa*kappa0; 

    g = kappa
end

% Chirped
if (strcmp(type,'chirped'))
    kappa0 = 100;
    temp = linspace(0,5,(Nz*stop)-(Nz*start)).'; % range from 0 to 1 over grating length
    kappa(round(Nz*start):round(Nz*start)+length(temp)-1) = temp;
    kappa = kappa0*kappa;
    g = kappa;
end

% Uniform Positive-Only Index Change
if (strcmp(type,'unipos'))
    kappa0 = 100;
    temp = linspace(0,3*pi,(Nz*stop)-(Nz*start)).'; % range from 0 to 1 over grating length
    kappa(round(Nz*start):round(Nz*start)+length(temp)-1) = abs(sin(temp));
    kappa = kappa0*kappa;
    g = kappa;
end


% Gaussian-Apodized Index Change
if (strcmp(type,'gaussian'))
    kappa0 = 100;
    temp = linspace(0,4*pi,(Nz*stop)-(Nz*start)).'; % range from 0 to 1 over grating length
    kappa(round(Nz*start):round(Nz*start)+length(temp)-1) = abs(sin(temp).*(4*exp(-((temp-6.5).^2)/2)));
    kappa = kappa0*kappa;
    g = kappa;
end

g = kappa;

