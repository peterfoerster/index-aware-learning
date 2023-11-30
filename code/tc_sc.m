pkg load symbolic

syms v_1 v_2 v_3 v_4 v_5 v_6 v_7 v_8 v_9 v_10 v_11 v_12 v_13 v_14 v_15 v_16
syms i_L_1 i_L_2 i_L_3 i_L_4 i_L_5 i_L_6 i_L_7
syms i_V_1 i_V_2

syms R_1 R_2 R_3 R_4 R_5 R_6 R_7
syms L_1 L_2 L_12 L_3 L_4 L_34 L_5 L_6 L_56 L_7
syms C_1 C_2 C_3 C_4 C_5 C_6 C_7 C_8 C_9
syms v_s_1 v_s_2

n_v   = 16;
n_R   = 7;
n_L   = 7;
n_C   = 9;
n_V   = 2;

A_R       = zeros(n_v,n_R);
A_R(1,1)  = 1;
A_R(2,1)  = -1;
A_R(8,2)  = 1;
A_R(9,2)  = -1;
A_R(3,3)  = 1;
A_R(4,3)  = -1;
A_R(10,4) = 1;
A_R(11,4) = -1;
A_R(5,5)  = 1;
A_R(6,5)  = -1;
A_R(12,6) = 1;
A_R(13,6) = -1;
A_R(12,7) = 1;
A_R(15,7) = -1;

A_L       = zeros(n_v,n_L);
A_L(2,1)  = 1;
A_L(3,1)  = -1;
A_L(9,2)  = 1;
A_L(10,2) = -1;
A_L(4,3)  = 1;
A_L(5,3)  = -1;
A_L(11,4) = 1;
A_L(12,4) = -1;
A_L(6,5)  = 1;
A_L(7,5)  = -1;
A_L(13,6) = 1;
A_L(14,6) = -1;
A_L(15,7) = 1;
A_L(16,7) = -1;

A_C       = zeros(n_v,n_C);
A_C(3,1)  = 1;
A_C(10,2) = 1;
A_C(3,2)  = -1;
A_C(5,3)  = 1;
A_C(12,4) = 1;
A_C(5,4)  = -1;
A_C(7,5)  = 1;
A_C(14,6) = 1;
A_C(7,6)  = -1;
A_C(10,7) = -1;
A_C(12,8) = -1;
A_C(16,9) = -1;

A_V      = zeros(n_v,n_V);
A_V(1,1) = 1;
A_V(8,2) = 1;

G      = sym(zeros(n_R,n_R));
G(1,1) = 1/R_1;
G(2,2) = 1/R_2;
G(3,3) = 1/R_3;
G(4,4) = 1/R_4;
G(5,5) = 1/R_5;
G(6,6) = 1/R_6;
G(7,7) = 1/R_7;

L      = sym(zeros(n_L,n_L));
L(1,1) = L_1;
L(1,2) = L_12;
L(2,1) = L_12;
L(2,2) = L_2;
L(3,3) = L_3;
L(3,4) = L_34;
L(4,3) = L_34;
L(4,4) = L_4;
L(5,5) = L_5;
L(5,6) = L_56;
L(6,5) = L_56;
L(6,6) = L_6;
L(7,7) = L_7;

C      = sym(zeros(n_C,n_C));
C(1,1) = C_1;
C(2,2) = C_2;
C(3,3) = C_3;
C(4,4) = C_4;
C(5,5) = C_5;
C(6,6) = C_6;
C(7,7) = C_7;
C(8,8) = C_8;
C(9,9) = C_9;

v_s = [v_s_1; v_s_2];

M = [[A_C*C*A_C.', zeros(n_v,n_L), zeros(n_v,n_V)]; [zeros(n_L,n_v), L, zeros(n_L,n_V)]; [zeros(n_V,n_v), zeros(n_V,n_L), zeros(n_V,n_V)];];
K = [[A_R*G*A_R.', A_L, A_V]; [-A_L.', zeros(n_L,n_L), zeros(n_L,n_V)]; [-A_V.', zeros(n_V,n_L), zeros(n_V,n_V)]];
f = [zeros(n_v,1); zeros(n_L,1); v_s];

% first basis functions
Q        = W        = null(M);
P        = V        = sym(zeros(n_v+n_L+n_V,14));
P(3,1)   = V(3,1)   = 1;
P(5,2)   = V(5,2)   = 1;
P(7,3)   = V(7,3)   = 1;
P(10,4)  = V(10,4)  = 1;
P(12,5)  = V(12,5)  = 1;
P(14,6)  = V(14,6)  = 1;
P(16,7)  = V(16,7)  = 1;
P(17,8)  = V(17,8)  = 1;
P(18,9)  = V(18,9)  = 1;
P(19,10) = V(19,10) = 1;
P(20,11) = V(20,11) = 1;
P(21,12) = V(21,12) = 1;
P(22,13) = V(22,13) = 1;
P(23,14) = V(23,14) = 1;

% first stage
Mt   = V.' * M * P;
Kt_P = V.' * K * P;
Kt_Q = V.' * K * Q;
ft   = V.' * f;
Kb_P = W.' * K * P;
Kb_Q = W.' * K * Q;
fb   = W.' * f;

% second basis functions
Qb = Wb = null(Kb_Q);
