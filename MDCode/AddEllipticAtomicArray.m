function AddEllipseAtomicArray(radX, radY, X0, Y0, VX0, VY0, InitDist, Temp, Type)
global C 
global x y AtomSpacing 
global nAtoms 
global AtomType Vx Vy Mass0 Mass1 Mass2

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

L = (2*radX - 1) * AtomSpacing;  
W = (2*radY - 1) * AtomSpacing;  

xp(1, :) = linspace(-L/2, L/2, 2*radX);
yp(1, :) = linspace(-W/2, W/2, 2*radY);

numAtoms = 0;
for i = 1:2*radX
    for j = 1:2*radY
        % Ellipse equation: (x/a)^2 + (y/b)^2 <= 1
        if (xp(i)/(radX*AtomSpacing))^2 + (yp(j)/(radY*AtomSpacing))^2 <= 1
            numAtoms = numAtoms+1;
            x(nAtoms + numAtoms) = xp(i);
            y(nAtoms + numAtoms) = yp(j);
        end
    end
end

x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + Y0;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);
    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;
end