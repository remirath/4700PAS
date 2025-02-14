% PA 4 Ryan Wood 
% 101195542

n_nodes = 50;                  % number of nodes on each side of square
size = 10;                      % size of square
n_iter = 100;                    % number of iterations

Vx = linspace(0,size,n_nodes);  % plate mesh x
Vy = linspace(0,size,n_nodes);  % plate mesh y

dx = size/(n_nodes-1);          % minimum x step
dy = size/(n_nodes-1);          % y step
n = 2:n_nodes-1;                % number of nodes on mesh in y
m = 2:n_nodes-1;                % in x

V = zeros(n_nodes, n_nodes);
E = zeros(n_nodes, n_nodes);

% M(n_iter) = struct('cdata',[],'colormap',[]);
% f = figure;
% f.Visible = 'off';
axis([0,size,0,size])                   % set axis to size
% ax.NextPlot = 'replaceChildren';

for j = 1:n_iter                        % update the equation for the number of iterations specified
    
    V(m,n) = (V(m+1,n)+V(m-1,n)+V(m,n+1)+V(m,n-1))/4;   % update equation for voltage
    E = -1*gradient(V,dx);                              % -1*gradient of voltage is electric field

    % Boundary conditions
    V(:,1) = 1.5;             % Left
    V(1,:) = 0;             % Bottom
    V(n_nodes,:) = 0;       % Top
    V(:,n_nodes) = 1.5;       % Right
    
    subplot(2,1,1)
    surf(Vx,Vy,V)           % surf plot the distribution of charge
    title('Voltage')
    view(0,90)              % set view to face-down
    shading interp
    colorbar
    % imboxfilt(V, 5, 'NormalizationFactor',1);
    drawnow

    subplot(2,1,2)
    surf(Vx,Vy,E)           % surface plot the electric field
    title('Electric Field')
    view(0,90)
    shading interp
    colorbar
    % imboxfilt(E, 5, 'NormalizationFactor',1);
    drawnow

    % M(j) = getframe;        % capture the plot as a movie frame
end

% f.Visible = 'on';
% movie(M,1,20);
