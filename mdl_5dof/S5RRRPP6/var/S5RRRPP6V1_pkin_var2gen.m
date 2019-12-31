% Umwandlung der Kinematikparameter von S5RRRPP6V1 zu S5RRRPP6
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRRPP6V1
%   pkin_var=[a4 a5 d1 theta4]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRPP6
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 theta4]
%
% Siehe auch: S5RRRPP6_structural_kinematic_parameters.m
function pkin_gen = S5RRRPP6V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([3, 4, 5, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(6) = 0.0; % d2
pkin_gen(7) = 0.0; % d3
