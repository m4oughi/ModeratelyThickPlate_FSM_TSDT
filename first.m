%{
Coder: Mohammad Foroughi(2024)

Method: Finite Strip Method
Subject: Buckling Analysis of Moderately Thick Isotropic Plate
Theory: Third Shear Deformation Theory

****************************************************************

Discriptions:
* Input file : 'inputs.xlsx'


*** Sections:
    A) inputs{
              E      : Modulus of Elasticity,
              noo    : Poisson's Ratio,
              G      : Shear's Modulus,
              total_a: Width of plate,
              b      : Length of plate,
              h      : Height of plate,
              m      : number of modes,
              n      : number of strips.
              ax, ay, axy: Directions of Applying Force
             }

    B) Requirements{
                    a: Width of Each Strip
                    D: Flextural Rigidity
                    Shape Functions --> Ni, Fi, Hi
                    }

    C) Calculating The Strain Matrix (Bm)
    D) Calculating Stiffness Matrix for Each Strip (k_local)
    E) Calculating Global Stiffness Matrix (k_global)
    F) Applying S-S Boundary Condition for Global Stiffness Matrix
    G) the displacements field (U, V, W)
    H) Calculating Geometric Stiffness Matrix (kg_local)
    I) Calculating Global Geometric Stiffness Matrix (kg_global)
    J) Applying S-S Boundary Condition for Global Geometric Stiffness Matrix
    K) Calculating Eigenvalue and Eigenvector
    L) Critical Buckling Load Parameter(minimum)
    M) Non-dimensionalization
%}

clc
clear all

%% A) Inputs
filename = 'inputs.xlsx';
inputs   = readtable(filename);

% --------  Modulus of Elasticity --------%
E = table2array(inputs(1, 3));

% --------     Poisson's Ratio    --------%
noo = table2array(inputs(2, 3));

% --------     Shear's Modulus    --------%
G = table2array(inputs(3, 3));

% --------        Dimensions      --------%
b       = table2array(inputs(4, 3));  % Length
total_a = table2array(inputs(5, 3));  % Width
h       = table2array(inputs(6, 3));  % height

% --------       Mode Number      --------%
m = table2array(inputs(7, 3));

% --------    Number of strips    --------%
n = table2array(inputs(8, 3));

% --------       Geometric        --------%
ax  = table2array(inputs(9, 3));
ay  = table2array(inputs(10, 3));
axy = table2array(inputs(11, 3));

kw_bar = table2array(inputs(12, 3));
ks_bar = table2array(inputs(13, 3));

%% B) Requirements
%
% --------   Width of Each Strip  --------%
a = total_a/n;

% --------   Flextural Rigidity  --------%
d = E/(1-noo^2);
D = d*[
    1     noo   0          0           0;
    noo   1     0          0           0;
    0     0    (1-noo)/2   0           0;
    0     0     0         (1-noo)/2    0;
    0     0     0          0          (1-noo)/2;
    ];

D0 = (E*h^3)/(12*(1-noo^2));

k1 = -4/(3*h^2);
k2 = -4/(h^2);

Kw = (kw_bar*D0)/(b^4);
Ks = (ks_bar*D0)/(b^2);

% ----  shape functions  ----%
syms x y z
N1 = ((2*x-a)*(x-a))/(a^2);
N2 = (4*x*(a-x))/(a^2);
N3 = (x*(2*x-a))/(a^2);

F1 = 1-(3*(x/a)^2)+(2*(x/a)^3);
F3 = (3*(x/a)^2)-(2*(x/a)^3);
H1 = x*(1-2*(x/a)+(x/a)^2);
H3 = x*(-(x/a)+(x/a)^2);
%}

%% C) Calculating Strain and displacement fields Matrices
%
B = []; N = []; U = []; V = []; W = [];
for i=1:1:m
    S = sin(i.*pi.*y./b);
    C = cos(i.*pi.*y./b);
    
    Bm(:, :, i)= [
        diff(N1, x)*S   0   (k1*z^3)*diff(F1, x, 2)*S   (k1*z^3)*diff(H1, x, 2)*S   (z+k1*z^3)*diff(N1, x)*S   0          diff(N2, x)*S   0   (z+k1*z^3)*diff(N2, x)*S   0          diff(N3, x)*S   0   (k1*z^3)*diff(F3, x, 2)*S   (k1*z^3)*diff(H3, x, 2)*S   (z+k1*z^3)*diff(N3, x)*S   0;
        0   N1*diff(C, y)   (k1*z^3)*F1*diff(S, y, 2)   (k1*z^3)*H1*diff(S, y, 2)   0   (z+k1*z^3)*N1*diff(C, y)          0   N2*diff(C, y)   0   (z+k1*z^3)*N2*diff(C, y)          0   N3*diff(C, y)   (k1*z^3)*F3*diff(S, y, 2)   (k1*z^3)*H3*diff(S, y, 2)   0   (z+k1*z^3)*N3*diff(C, y);
        N1*diff(S, y)   diff(N1, x)*C   2*(k1*z^3)*diff(F1, x)*diff(S, y)   2*(k1*z^3)*diff(H1, x)*diff(S, y)   (z+k1*z^3)*N1*diff(S, y)   (z+k1*z^3)*diff(N1, x)*C          N2*diff(S, y)   diff(N2, x)*C   (z+k1*z^3)*N2*diff(S, y)   (z+k1*z^3)*diff(N2, x)*C          N3*diff(S, y)   diff(N3, x)*C   2*(k1*z^3)*diff(F3, x)*diff(S, y)   2*(k1*z^3)*diff(H3, x)*diff(S, y)   (z+k1*z^3)*N3*diff(S, y)   (z+k1*z^3)*diff(N3, x)*C;
        0   0   (1+k2*z^2)*F1*diff(S, y)   (1+k2*z^2)*H1*diff(S, y)   0   (1+k2*z^2)*N1*C          0   0   0   (1+k2*z^2)*N2*C          0   0   (1+k2*z^2)*F3*diff(S, y)   (1+k2*z^2)*H3*diff(S, y)   0   (1+k2*z^2)*N3*C;
        0   0   (1+k2*z^2)*diff(F1, x)*S   (1+k2*z^2)*diff(H1, x)*S   (1+k2*z^2)*N1*S   0          0   0   (1+k2*z^2)*N2*S   0          0   0   (1+k2*z^2)*diff(F3, x)*S   (1+k2*z^2)*diff(H3, x)*S   (1+k2*z^2)*N3*S   0;
        ];
    
    Nm(:, :, i) = [
        N1*S   0   (k1*z^3)*diff(F1, x)*S   (k1*z^3)*diff(H1, x)*S   (z+k1*z^3)*N1*S   0          N2*S   0   (z+k1*z^3)*N2*S   0          N3*S   0   (k1*z^3)*diff(F3, x)*S   (k1*z^3)*diff(H3, x)*S   (z+k1*z^3)*N3*S   0;
        0   N1*C   (k1*z^3)*F1*diff(S, y)   (k1*z^3)*H1*diff(S, y)   0   (z+k1*z^3)*N1*C          0   N2*C   0   (z+k1*z^3)*N2*C          0   N3*C   (k1*z^3)*F3*diff(S, y)   (k1*z^3)*H3*diff(S, y)   0   (z+k1*z^3)*N3*C;
        0   0   F1*S   H1*S   0   0   0   0   0   0   0   0   F3*S   H3*S   0   0;
        ];
    
    Um(:, :, i) = [N1*S 0 (k1*z^3)*diff(F1, x)*S (k1*z^3)*diff(H1, x)*S (z+k1*z^3)*N1*S 0 N2*S 0 (z+k1*z^3)*N2*S 0 N3*S 0 (k1*z^3)*diff(F3, x)*S (k1*z^3)*diff(H3, x)*S (z+k1*z^3)*N3*S 0];
    Vm(:, :, i) = [0 N1*C (k1*z^3)*F1*diff(S, y) (k1*z^3)*H1*diff(S, y) 0 (z+k1*z^3)*N1*C 0 N2*C 0 (z+k1*z^3)*N2*C 0 N3*C (k1*z^3)*F3*diff(S, y) (k1*z^3)*H3*diff(S, y) 0 (z+k1*z^3)*N3*C];
    Wm(:, :, i) = [0 0 F1*S H1*S 0 0 0 0 0 0 0 0 F3*S H3*S 0 0];
    
    B = [B, Bm(:, :, i)];
    N = [N, Nm(:, :, i)];
    U = [U, Um(:, :, i)];
    V = [V, Vm(:, :, i)];
    W = [W, Wm(:, :, i)];
end
disp("B, N, U, V, W is generated!");
%}

%%  D) Calculating Stiffness Matrix for Each Strip (16*16)
%
k__local = int(int(int( B.' * D * B, x, 0, a), y, 0, b), z, -0.5*h, 0.5*h);
k_local = double(k__local);
disp("Local stiffness matrix is calculated!")
%}

%% Calculating Elastic foundation impact
%
% Winkler
kw_local = int(int( subs((N.'*Kw*N), z, 0), x, 0, a), y, 0, b);
kw_local = double(kw_local);
disp("Winkler stiffness matrix is calculated!")
%}

%
% Pasternak
Nx = diff(N, x);
Ny = diff(N, y);
ksx_local = int(int( subs(Nx.'*Ks*Nx, z, 0), x, 0, a), y, 0, b);
ksy_local = int(int( subs(Ny.'*Ks*Ny, z, 0), x, 0, a), y, 0, b);
ks_local = double(ksx_local + ksy_local);
disp("Pasternak stiffness matrix is calculated!")
%}

%
kLocal = k_local + kw_local + ks_local;
disp("Sum of local, winkler, and pasternak stiffness matrices is calculated!")
%}

%% E) Calculating Global Stiffness Matrix
%
kltokg1 = zeros(16, 16, m, m);
for i = 1:1:m
    for j = 1:1:m
        kltokg1(:, :, i, j) = kLocal(16*i-15:16*i, 16*j-15:16*j);
    end
end

kltokg2 = zeros(10*n+6, 10*n+6, m, m);
for i = 1:1:m
    for j = 1:1:m
        for k = 1:1:n
            kltokg3 = zeros(10*n+6, 10*n+6, m, m);
            kltokg3(10*k-9:10*k+6, 10*k-9:10*k+6, i, j) = kltokg1(:, :, i, j);
            kltokg2(:, :, i, j) = kltokg2(:, :, i, j) + kltokg3(:, :, i, j);
        end
    end
end

k_global = zeros(m*(10*n+6), m*(10*n+6));
for i = 1:1:m
    for j = 1:1:m
        k_global((10*n+6)*i-(10*n+5):(10*n+6)*i, (10*n+6)*j-(10*n+5):(10*n+6)*j) = kltokg2(:, :, i, j);
    end
end
disp("Global stiffness matrix is calculated!")
%}


%% G) the displacements field (16*1)
%
U = U.';
V = V.';
W = W.';
%}

%% H) Calculating Geometric Stiffness Matrix (16*16)
%
kg__local = int( int( int( ax*( diff(U,x).*diff(U,x).' + diff(V,x).*diff(V,x).' + diff(W,x).*diff(W,x).' ) + ay*( diff(U, y).* diff(U,y).' + diff(V,y).* diff(V,y).' + diff(W,y).* diff(W,y).' ) + axy*( diff(U,x).*diff(U,y).' + diff(U,x).*diff(U,y).' + diff(V,x).*diff(V,y).' + diff(V,x).*diff(V,y).' + diff(W,x).* diff(W,y).' + diff(W,x).*diff(W,y).' ), x, 0, a), y, 0, b), z, -0.5*h, 0.5*h);
kg_local = double(kg__local);
disp("Gemoetry stiffness matrix is calculated!")
%}

%% I) Calculating Global Stiffness Matrix (10*n+6)*(10*n+6)
%
kgtokgg1 = zeros(16, 16, m, m);
for i = 1:1:m
    for j = 1:1:m
        kgtokgg1(:, :, i, j) = kg_local(16*i-15:16*i, 16*j-15:16*j);
    end
end

kgtokgg2 = zeros(10*n+6, 10*n+6, m, m);
for i = 1:1:m
    for j = 1:1:m
        for k = 1:1:n
            kgtokgg3 = zeros(10*n+6, 10*n+6, m, m);
            kgtokgg3(10*k-9:10*k+6, 10*k-9:10*k+6, i, j) = kgtokgg1(:, :, i, j);
            kgtokgg2(:, :, i, j) = kgtokgg2(:, :, i, j) + kgtokgg3(:, :, i, j);
        end
    end
end

kg_global = zeros(m*(10*n+6), m*(10*n+6));
for i = 1:1:m
    for j = 1:1:m
        kg_global((10*n+6)*i-(10*n+5):(10*n+6)*i, (10*n+6)*j-(10*n+5):(10*n+6)*j) = kgtokgg2(:, :, i, j);
    end
end
disp("Global geometry stiffness matrix is calculated!")
%}

%% F) Applying S-S Boundary Condition for Global Stiffness Matrix (10*n)*(10*n)
rm = [];
%
%S-S
for i = m:-1:1
    rm = [rm, (10*n+6)*i];
    rm = [rm, (10*n+6)*i-3];
    rm = [rm, (10*n+6)*i-4];
    %rm = [rm, (10*n+6)*i-5];
    rm = [rm, (10*n+6)*(i-1)+6];
    rm = [rm, (10*n+6)*(i-1)+3];
    rm = [rm, (10*n+6)*(i-1)+2];
    %rm = [rm, (10*n+6)*(i-1)+1];
end
%}
%{
% C-C
for i = m:-1:1
    rm = [rm, (10*n+6)*i];
    rm = [rm, (10*n+6)*i-1];
    rm = [rm, (10*n+6)*i-2];
    rm = [rm, (10*n+6)*i-3];
    rm = [rm, (10*n+6)*i-4];
    rm = [rm, (10*n+6)*i-5];
    rm = [rm, (10*n+6)*(i-1)+6];
    rm = [rm, (10*n+6)*(i-1)+5];
    rm = [rm, (10*n+6)*(i-1)+4];
    rm = [rm, (10*n+6)*(i-1)+3];
    rm = [rm, (10*n+6)*(i-1)+2];
    rm = [rm, (10*n+6)*(i-1)+1];
end

%{
% C-F
for i = m:-1:1
    rm = [rm, (10*n+6)*i];
    rm = [rm, (10*n+6)*i-1];
    rm = [rm, (10*n+6)*i-2];
    rm = [rm, (10*n+6)*i-3];
    rm = [rm, (10*n+6)*i-4];
    rm = [rm, (10*n+6)*i-5];
end
%}

%
%S-F
for i = m:-1:1
    rm = [rm, (10*n+6)*i];
    rm = [rm, (10*n+6)*i-3];
    rm = [rm, (10*n+6)*i-4];
    rm = [rm, (10*n+6)*i-5];
end
%}
%{
% S-C
for i = m:-1:1
    rm = [rm, (10*n+6)*i];
    rm = [rm, (10*n+6)*i-3];
    rm = [rm, (10*n+6)*i-4];
    rm = [rm, (10*n+6)*i-5];
    rm = [rm, (10*n+6)*(i-1)+6];
    rm = [rm, (10*n+6)*(i-1)+5];
    rm = [rm, (10*n+6)*(i-1)+4];
    rm = [rm, (10*n+6)*(i-1)+3];
    rm = [rm, (10*n+6)*(i-1)+2];
    rm = [rm, (10*n+6)*(i-1)+1];
end
%}
s_rm = size(rm);
for i = 1:1:s_rm(2)
    k_global(rm(i), :) = [];
    kg_global(rm(i), :) = [];
    
    k_global(:, rm(i)) = [];
    kg_global(:, rm(i)) = [];
end
%}


%% K) Calculating Eigenvalue and Eigenvector
[Vec, Landa] = eig(k_global, kg_global);

%% L) Critical Buckling Load Parameter(minimum)
Max_D = max(Landa);
dddd = min(Max_D);


%% M) Non-dimensionalization
dddd*((total_a)^2)*12*(1-noo^2)*(b^2) / ((pi^2)*(h^2)*E)
dddd*((total_a)^2)*12*(1-noo^2)*(b^2) / ((h^2)*E)

%}
