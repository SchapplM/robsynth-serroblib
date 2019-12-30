% Umwandlung der Kinematikparameter von S5PRRPR7V2 zu S5PRRPR7
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S5PRRPR7V2
%   pkin_var=[a2 a3 a4 a5 alpha2 d2 d3 d5 theta1]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5PRRPR7
%   pkin_gen=[a2 a3 a4 a5 alpha2 d2 d3 d5 theta1 theta4]
%
% Siehe auch: S5PRRPR7_structural_kinematic_parameters.m
function pkin_gen = S5PRRPR7V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 9]) = pkin_var;

pkin_gen(10) = pi/2; % theta4
