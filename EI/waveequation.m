% PA 5 Ryan Wood 101195542
% clear all
set(0, 'DefaultFigureWindowStyle','docked')

nx = 50;
ny = 50;                  % number of nodes on each side of square
V = zeros(nx,ny);
G = sparse(nx*ny,nx*ny);
Inclusion = 0;

for i = 1:nx
    for j = 1:ny

        n = j + (i-1)*ny
        nxm = j + (i-2)*ny
        nxp = j +i*ny
        nym = (j-1) + (i-1)*ny
        nyp = (j+1) + (i-1)*ny

        % V(n) = (1/4)*(V(nxp,n)+V(nxm,n)+V(n,nyp)+V(n,nym));    %update equation

        if(i == 1 || i == nx || j == 1 || j == ny)
            G(n,n) = 1;
        elseif(i > 10 && i < 20 && j > 10 && j < 20)
            G(n,n) = -2;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        else 
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end

    end
end

figure('Name','Matrix')
spy(G)

nmodes = 20;
[E,D] = eigs(G,nmodes,'SM');

figure('Name','Eigenvalues')
plot(diag(D),'*');

np = ceil(sqrt(nmodes))
figure('name','Modes')
for k = 1:nmodes
    M = E(:,k);
    for i = 1:nx
        for j = 1:ny
            n = i + (j-1)*nx;
            V(i,j) = M(n);
        end
        subplot(np,np,k), surf(V,'linestyle','none')
        title(['EV = ' num2str(D(k,k))])
    end
end

