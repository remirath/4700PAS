% This example shows how to calculate and plot both the
% fundamental TE and TM eigenmodes of an example 3-layer ridge
% waveguide using the full-vector eigenmode solver.  

% Refractive indices:
n1 = 3.34;          % Lower cladding
n2 = 3.44;          % Core
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding

% Horizontal dimensions:
rh = 1.1;           % Ridge height
rw = 1.0;           % Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 5;         % number of modes to compute
te = 1;             % To activate TE mode set to 1, otherwise 0
tm = 1;             % To activate TM mode set to 1, otherwise 0

[x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,n2,n3],[h1,h2,h3], ...
                                            rh,rw,side,dx,dy); 
i = 1;
j = 1;

% Consider both the fundamental TM and TE modes

[Hx,Hy,neff] = wgmodes(lambda,n2,nmodes,dx,dy,eps,'000A'); % TE mode
[Hx2,Hy2,neff2] = wgmodes(lambda,n2,nmodes,dx,dy,eps,'000S'); % TM mode
while i <= nmodes

% TE mode
if te==1
fprintf(1,'neff = %.6f\n',neff);

figure(j);
subplot(121);
contourmode(x,y,Hx(:,:,i));
title('Hx (TE mode)'); xlabel('x'); ylabel('y'); 
for v = edges, line(v{:}); end

subplot(122);
contourmode(x,y,Hy(:,:,i));
title('Hy (TE mode)'); xlabel('x'); ylabel('y'); 
for v = edges, line(v{:}); end
end
% TM mode
if tm==1
fprintf(1,'neff = %.6f\n',neff2);

figure(j+1);
subplot(121);
contourmode(x,y,Hx2(:,:,i));
title('Hx (TM mode)'); xlabel('x'); ylabel('y'); 
for v = edges, line(v{:}); end

subplot(122);
contourmode(x,y,Hy2(:,:,i));
title('Hy (TM mode)'); xlabel('x'); ylabel('y'); 
for v = edges, line(v{:}); end
end

i = i + 1;
j = j + 2;

end

% iterate ridge half-width from 0.325 to 1.0

% n = 1;
% n2 = 0.325;
% for n = 1:11
%     [x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,n2,n3],[h1,h2,h3], ...
%                                             rh,rw,side,dx,dy);
%     [Hx,Hy,neff] = wgmodes(lambda,n2,nmodes,dx,dy,eps,'000A');
%     fprintf(1,'neff = %.6f\n',neff);
%     fprintf(1,'rw = %.6f\n',rw);
% 
%     figure(n);
%     subplot(121);
%     contourmode(x,y,Hx);
%     title('Hx (TE mode)'); xlabel('x'); ylabel('y'); 
%     for v = edges, line(v{:}); end
% 
%     subplot(122);
%     contourmode(x,y,Hy);
%     title('Hy (TE mode)'); xlabel('x'); ylabel('y'); 
%     for v = edges, line(v{:}); end
%     n2 = n2+0.0.0675;
%     n = n + 1;
% 
% end

% iterate ridge index from 3.305 to 3.44

% n = 1;
% n2 = 3.305;
% for n = 1:11
%     [x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,n2,n3],[h1,h2,h3], ...
%                                             rh,rw,side,dx,dy);
%     [Hx,Hy,neff] = wgmodes(lambda,n2,nmodes,dx,dy,eps,'000A');
%     fprintf(1,'neff = %.6f\n',neff);
%     fprintf(1,'n2 = %.6f\n',n2);
% 
%     figure(n);
%     subplot(121);
%     contourmode(x,y,Hx);
%     title('Hx (TE mode)'); xlabel('x'); ylabel('y'); 
%     for v = edges, line(v{:}); end
% 
%     subplot(122);
%     contourmode(x,y,Hy);
%     title('Hy (TE mode)'); xlabel('x'); ylabel('y'); 
%     for v = edges, line(v{:}); end
%     n2 = n2 + 0.0135;
%     n = n + 1;
% 
% end