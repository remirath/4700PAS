function g = grating(start,stop,type,kappa0)

% This function models an assortment of Fiber Bragg Grating (FBG)
% structures which can be used in the simulation of a 1D optical waveguide.
% Returns a Nx1 array which represents the coupling coefficient between
% waveguide indexes at each point. The overall coupling coefficients along
% the waveguide will determine the filter response and center frequency for
% reflection (Bragg wavelength). 
% 
% Inputs
% 
% start, stop - The grating start and stop points as a fraction of the total waveguide
% length from 0 to 1.
%
% type - The type of grating desired, uniform, chirped, unipos, and
% gaussian. Leave blank for no grating.
%
% kappa0 - The characteristic coupling coefficient. Can be altered to
% change the Bragg frequency and bandwidth of each grating.
%

L = 1000e-6*1e2;             % cm length of waveguide
Nz = 500;                    % number of points along z axis
z = linspace(0,L,Nz).';      % Nz points, Nz-1 segments
kappa = zeros(size(z));      % kappa as matrix

% Uniform
if (strcmp(type,'uniform'))
    temp = linspace(1,1,(Nz*stop)-(Nz*start)).';
    %
    % y = 1
    kappa(round(Nz*start):round(Nz*start)+length(temp)-1) = temp;
    kappa = kappa*kappa0; 
end

% Chirped
if (strcmp(type,'chirped'))
    temp = linspace(0,1,(Nz*stop)-(Nz*start)).'; % range from 0 to 1 over grating length
    temp_sin = linspace(0,2*pi,(Nz*stop)-(Nz*start)).';

    % The chirped grating is a uniform positive-only grating scaled
    % linearly from 0 to 1. This results in large positive values followed
    % by zero in a periodic fashion.
    % y = x*abs(sin(2x))
    kappa(round(Nz*start):round(Nz*start)+length(temp)-1) = temp.*abs(sin(temp_sin.*2));
    kappa = kappa0*kappa;
end

% Uniform Positive-Only Index Change
if (strcmp(type,'unipos'))
    temp = linspace(0,2*pi,(Nz*stop)-(Nz*start)).'; % range from 0 to 1 over grating length

    % The uniform positive-only grating is the absolute value of a sin
    % function, ranging from large positive values to 0 periodically.
    % y = abs(sin(3x))
    kappa(round(Nz*start):round(Nz*start)+length(temp)-1) = abs(sin(temp.*3));
    kappa = kappa0*kappa;
end


% Gaussian-Apodized Index Change
if (strcmp(type,'gaussian'))
    temp = linspace(0,2*pi,(Nz*stop)-(Nz*start)).'; % range from 0 to 4pi over grating length

    % The gaussian-apodized function shape is a gaussian pulse modulated by
    % a sin function.
    % y = abs(sin(3x)*4e^-(x-6.5)^2/2)
    kappa(round(Nz*start):round(Nz*start)+length(temp)-1) = abs(sin(3.*temp).*(4*exp(-((temp-6.5).^2)/2)));
    kappa = kappa0*kappa;
end

g = kappa;      % return the array. array is 0 if grating is left unspecified

