% Umwandlung der Kinematikparameter von S5PRRPR7V4 zu S5PRRPR7
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5PRRPR7V4
%   pkin_var=[a2 a4 a5 d2 d5 theta1]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5PRRPR7
%   pkin_gen=[a2 a3 a4 a5 alpha2 d2 d3 d5 theta1 theta4]
%
% Siehe auch: S5PRRPR7_structural_kinematic_parameters.m
function pkin_gen = S5PRRPR7V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 3, 4, 6, 8, 9]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(5) = pi/2; % alpha2
pkin_gen(7) = 0.0; % d3
pkin_gen(10) = 0.0; % theta4
