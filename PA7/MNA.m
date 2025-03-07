G1 = 1;
G2 = 0.5;
G3 = 0.1;
G4 = 10;
G0 = 0.001;
C = 0.25;
L = 0.2;
a = 100;
vin = 1;
v1 = 1;
v2 = 1;
v3 = 1;
va = 1;
v4 = 1;
v0 = 1;
iL = 1;
iA = 1;
iS = 1;


G = [ 1      0      0     0     0    0      0    0   0   
      G1    -G1     0     0     0    0      1    0   0   
     -G1   G1+G2    0     0     0    0      0   -1   0   
      0      -1    1-G3   0     0    0      0    1   0   
      0       0     0    a*iL  -1    0      0    0   0      
      0       0     0     0    -G4   G4     0    0   1   
      0       0     0     0    -G4  G4-G0   0    0   0   
      0       1     0     0     0    0      0    0   0   
      0       0     -1    1     0    0      0    0   0   
      0       0     0     0     0    1      0    0   0     
    ];

C = [ 0      0      0     0     0    0      0    0   0   
      C      -C      0     0     0    0      0    0   0   
     -C      C      0     0     0    0      0    0   0   
      0       0     0     0     0    0      0    0   0   
      0       0     0     0     0    0      0    0   0   
      0       0     0     0     0    0      0    0   0   
      0       0     0     0     0    0      0    0   0   
      0       0     0     0     0    0      0    0   0   
      0       0     0     0     0    0      0    -L   0   
      0       0     0     0     0    0      0    0   0   
    ]; 

V = [ 
      v1
      v2
      v3
      va
      v4 
      v0
      iS  
      iL
      iA 
    ];

F = [1
     0
     0
     0
     0
     0
     0
     0
     0
     0
    ];

    figure('name','Fields')             % set up plotting
    plot(0:20,v0,'r');hold on 
    plot(0:20,v3,'b');          % plot real part of forward envelope
    xlabel('Vin (V)')
    ylabel('V0')
    hold off


for i = -10:10
        
        F = C+(G*V);

        plot(0:20,v0,'r'); hold on
        plot(0:20,v3,'b');
        hold off
end