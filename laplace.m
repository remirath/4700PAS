% PA 4 Ryan Wood 

n_nodes = 100;                  % number of nodes on each side of square
size = 10;                      % size of square
n_iter = 100;                    % number of iterations

Vx = linspace(0,size,n_nodes);  % 
Vy = linspace(0,size,n_nodes);

dx = size/(n_nodes-1);
dy = size/(n_nodes-1);
n = 2:n_nodes-1;
m = 2:n_nodes-1;

V = zeros(n_nodes, n_nodes);

M(n_iter) = struct('cdata',[],'colormap',[]);
f = figure;
% f.Visible = 'off';
axis manual
axis([0,10,0,10])
hold on

for j = 1:n_iter
    V(m,n) = ((V(m+1,n)-2*V(m,n)+V(m-1,n))/dx^2)+((V(m,n-1)-2*V(m,n)+V(m,n+1))/dy^2)

    V(:,1) = 1;
    V(1,:) = 1;
    V(n_nodes,:) = 1;
    V(:,n_nodes) = 1;

    surf(Vx,Vy,V)

    M(j) = getframe;
end

f.Visible = 'on';

