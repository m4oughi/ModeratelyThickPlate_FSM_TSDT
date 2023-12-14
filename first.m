%{
Coder: Mohammad Foroughi(Jun 2020)

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
clear

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

% --------    Number of Nodal lines    --------%
n = table2array(inputs(8, 3));

% --------       Geometric        --------%
ax  = table2array(inputs(9, 3));
ay  = table2array(inputs(10, 3));
axy = table2array(inputs(11, 3));



%% B) Requirements
% --------   Width of Each Strip  --------%
a = total_a/n;

% --------   Flextural Rigidity  --------%
% COMPLETELY CHECKED
d = E/(1-noo^2);
D = d*[
    1     noo   0          0           0;
    noo   1     0          0           0;
    0     0    (1-noo)/2   0           0;
    0     0     0         (1-noo)/2    0;
    0     0     0          0          (1-noo)/2;
    ];


% ----  shape functions  ----%
syms z x y
N1 = (2*x-a)*(x-a)/(a^2);
N2 = (4*x*(a-x)/(a^2));
N3 = (x*(2*x-a))/(a^2);

F1 = 1-(3*(x/a)^2)+(2*(x/a)^3);
F3 = (3*(x/a)^2)-(2*(x/a)^3);
H1 = x*(1-2*(x/a)+(x/a)^2);
H3 = x*(-(x/a)+(x/a)^2);

S = sym(zeros(1, m));
C = sym(zeros(1, m));
for i=1:1:m
    S(1, i) = sin(i*pi*y/b);
    C(1, i) = cos(i*pi*y/b);
end


%% C) Calculating Strain Matrix (5*16)
k1 = -4/(3*h^2);
k2 = -4/(h^2);

Bm = sym(zeros(5, 16, m));
for i=1:1:m
    Bm(:, :, i) = [
        diff(N1, x)*S(1, i)   0   (k1*z^3)*diff(F1, x, 2)*S(1, i)   (k1*z^3)*diff(H1, x, 2)*S(1, i)   (z+k1*z^3)*diff(N1, x)*S(1, i)   0          diff(N2, x)*S(1, i)   0   (z+k1*z^3)*diff(N2, x)*S(1, i)   0          diff(N3, x)*S(1, i)   0   (k1*z^3)*diff(F3, x, 2)*S(1, i)   (k1*z^3)*diff(H3, x, 2)*S(1, i)   (z+k1*z^3)*diff(N3, x)*S(1, i)   0;
        0   N1*diff(C(1, i), y)   (k1*z^3)*F1*diff(S(1, i), y, 2)   (k1*z^3)*H1*diff(S(1, i), y, 2)   0   (z+k1*z^3)*N1*diff(C(1, i), y)          0   N2*diff(C(1, i), y)   0   (z+k1*z^3)*N2*diff(C(1, i), y)          0   N3*diff(C(1, i), y)   (k1*z^3)*F3*diff(S(1, i), y, 2)   (k1*z^3)*H3*diff(S(1, i), y, 2)   0   (z+k1*z^3)*N3*diff(C(1, i), y);
        N1*diff(S(1, i), y)   diff(N1, x)*C(1, i)   2*(k1*z^3)*diff(F1, x)*diff(S(1, i), y)   2*(k1*z^3)*diff(H1, x)*diff(S(1, i), y)   (z+k1*z^3)*N1*diff(S(1, i), y)   (z+k1*z^3)*diff(N1, x)*C(1, i)          N2*diff(S(1, i), y)   diff(N2, x)*C(1, i)   (z+k1*z^3)*N2*diff(S(1, i), y)   (z+k1*z^3)*diff(N2, x)*C(1, i)          N3*diff(S(1, i), y)   diff(N3, x)*C(1, i)   2*(k1*z^3)*diff(F3, x)*diff(S(1, i), y)   2*(k1*z^3)*diff(H3, x)*diff(S(1, i), y)   (z+k1*z^3)*N3*diff(S(1, i), y)   (z+k1*z^3)*diff(N3, x)*C(1, i);
        0   0   (1+k2*z^2)*F1*diff(S(1, i), y)   (1+k2*z^2)*H1*diff(S(1, i), y)   0   (1+k2*z^2)*N1*C(1, i)          0   0   0   (1+k2*z^2)*N2*C(1, i)          0   0   (1+k2*z^2)*F3*diff(S(1, i), y)   (1+k2*z^2)*H3*diff(S(1, i), y)   0   (1+k2*z^2)*N3*C(1, i);
        0   0   (1+k2*z^2)*diff(F1, x)*S(1, i)   (1+k2*z^2)*diff(H1, x)*S(1, i)   (1+k2*z^2)*N1*S(1, i)   0          0   0   (1+k2*z^2)*N2*S(1, i)   0          0   0   (1+k2*z^2)*diff(F3, x)*S(1, i)   (1+k2*z^2)*diff(H3, x)*S(1, i)   (1+k2*z^2)*N3*S(1, i)   0;
        ];
end


%%  D) Calculating Stiffness Matrix for Each Strip (16*16)

k__local = sym(zeros(16, 16, n, m));
for i = 1:1:m
    for j = 1:1:n
        k__local(:, :, j, i) = int(int(int( Bm(:, :, i)' * D * Bm(:, :, i), x, 0, a), y, 0, b), z, -0.5*h, 0.5*h);
    end
end
k_local = double(k__local);



%% E) Calculating Global Stiffness Matrix (10*n+6)*(10*n+6)
k__global = sym(zeros(10*n+6, 10*n+6, m));
for i=1:1:m
    for j = 1:1:n
        k___global = sym(zeros(10*n+6, 10*n+6));
        k___global(10*j-9:10*j+6, 10*j-9:10*j+6) = k_local(:, :, j, i);
        k__global(:, :, i) = k__global(:, :, i) + k___global;
    end
end
k_global = double(k__global);



%% F) Applying S-S Boundary Condition for Global Stiffness Matrix (10*n)*(10*n)
% END
k_global(10*n+6, :, :) = [];    % remove W from last strip(row)
k_global(:, 10*n+6, :) = [];    % remove W from last strip(column)
k_global(10*n+3, :, :) = [];    % remove V from last strip(row)
k_global(:, 10*n+3, :) = [];    % remove V from last strip(column)
k_global(10*n+2, :, :) = [];    % remove U from last strip(row)
k_global(:, 10*n+2, :) = [];    % remove U from last strip(column)

% FIRST
k_global(6, :, :) = [];        % remove W from first strip(row)
k_global(:, 6, :) = [];        % remove W from first strip(column)
k_global(3, :, :) = [];        % remove V from first strip(row)
k_global(:, 3, :) = [];        % remove V from first strip(column)
k_global(2, :, :) = [];        % remove U from first strip(row)
k_global(:, 2, :) = [];        % remove U from first strip(column)



%% G) the displacements field (16*1)
U = sym(zeros(16, 1, m));
V = sym(zeros(16, 1, m));
W = sym(zeros(16, 1, m));
for i=1:1:m
    U(:, 1, i) = [N1*S(1, i);   0;   (k1*z^3)*diff(F1, x)*S(1, i);   (k1*z^3)*diff(H1, x)*S(1, i);   (z+k1*z^3)*N1*S(1, i);   0;          N2*S(1, i);   0;   (z+k1*z^3)*N2*S(1, i);   0;          N3*S(1, i);   0;   (k1*z^3)*diff(F3, x)*S(1, i);   (k1*z^3)*diff(H3, x)*S(1, i);   (z+k1*z^3)*N3*S(1, i);   0];
    V(:, 1, i) = [0;   N1*C(1, i);   (k1*z^3)*F1*diff(S(1, i), y);   (k1*z^3)*H1*diff(S(1, i), y);   0;   (z+k1*z^3)*N1*C(1, i);          0;   N2*C(1, i);   0;   (z+k1*z^3)*N2*C(1, i);          0;   N3*C(1, i);   (k1*z^3)*F3*diff(S(1, i), y);   (k1*z^3)*H3*diff(S(1, i), y);   0;   (z+k1*z^3)*N3*C(1, i)];
    W(:, 1, i) = [0;   0;   F1*S(1, i);   H1*S(1, i);   0;   0;   0;   0;   0;   0;   0;   0;   F3*S(1, i);   H3*S(1, i);   0;   0];
end

%% H) Calculating Geometric Stiffness Matrix (16*16)
kg__local = sym(zeros(16, 16, n, m));
for i=1:1:m
    for j = 1:1:n
        kg__local(:, :, j, i) = int( int( int( ax*( diff(U(:, :, i),x).*diff(U(:, :, i),x)' + diff(V(:, :, i),x).*diff(V(:, :, i),x)' + diff(W(:, :, i),x).*diff(W(:, :, i),x)' ) + ay*( diff(U(:, :, i), y).* diff(U(:, :, i),y)' + diff(V(:, :, i),y).* diff(V(:, :, i),y)' + diff(W(:, :, i),y).* diff(W(:, :, i),y)' ) + axy*( diff(U(:, :, i),x).*diff(U(:, :, i),y)' + diff(U(:, :, i),x).*diff(U(:, :, i),y)' + diff(V(:, :, i),x).*diff(V(:, :, i),y)' + diff(V(:, :, i),x).*diff(V(:, :, i),y)' + diff(W(:, :, i),x).* diff(W(:, :, i),y)' + diff(W(:, :, i),x).*diff(W(:, :, i),y)' ), x, 0, a), y, 0, b), z, -0.5*h, 0.5*h);
    end
end
kg_local = double(kg__local);



%% I) Calculating Global Stiffness Matrix (10*n+6)*(10*n+6)
kg__global = sym(zeros(10*n+6, 10*n+6, m));
for i=1:1:m
    for j = 1:1:n
        kg___global = sym(zeros(10*n+6, 10*n+6));
        kg___global(10*j-9:10*j+6, 10*j-9:10*j+6) = kg__local(:, :, j, i);
        kg__global(:, :, i) = kg__global(:, :, i) + kg___global;
    end
end
kg_global = double(kg__global);



%% J) S-S Boundary Condition for Global Geometric Stiffness Matrix (10*n)*(9*n)

% END
kg_global(10*n+6, :, :) = [];   % remove W from last strip(row)
kg_global(:, 10*n+6, :) = [];   % remove W from last strip(column)
kg_global(10*n+3, :, :) = [];   % remove V from last strip(row)
kg_global(:, 10*n+3, :) = [];   % remove V from last strip(column)
kg_global(10*n+2, :, :) = [];   % remove U from last strip(row)
kg_global(:, 10*n+2, :) = [];   % remove U from last strip(column)

% FIRST
kg_global(6, :, :) = [];       % remove W from first strip(row)
kg_global(:, 6, :) = [];       % remove W from first strip(column)
kg_global(3, :, :) = [];       % remove V from first strip(row)
kg_global(:, 3, :) = [];       % remove V from first strip(column)
kg_global(2, :, :) = [];       % remove U from first strip(row)
kg_global(:, 2, :) = [];       % remove U from first strip(column)



%% K) Calculating Eigenvalue and Eigenvector
Vec = zeros(10*n , 10*n, m);
Landa = zeros(10*n, 10*n, m);
for i=1:1:m
    [Vec(:, :, i), Landa(:, :, i)] = eig(k_global(:, :, i), kg_global(:, :, i));
end



%% L) Critical Buckling Load Parameter(minimum)
Max_D = zeros(1, 10*n, m);
dddd = zeros(1, m);
for i=1:1:m
    Max_D(:, :, i) = max(Landa(:, :, i));
    dddd(:, i) = min(Max_D(:, :, i));
end



%% M) Non-dimensionalization
for i=1:1:m
    dddd(1, i)*((total_a)^2)*12*(1-noo^2) / ((h^2)*E)
end
%}
