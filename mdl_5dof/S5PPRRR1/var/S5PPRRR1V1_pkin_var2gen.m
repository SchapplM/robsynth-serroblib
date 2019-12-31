% Umwandlung der Kinematikparameter von S5PPRRR1V1 zu S5PPRRR1
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PPRRR1V1
%   pkin_var=[a2 a3 a4 d3 d4 theta1 theta2]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PPRRR1
%   pkin_gen=[a2 a3 a4 a5 d3 d4 d5 theta1 theta2]
%
% Siehe auch: S5PPRRR1_structural_kinematic_parameters.m
function pkin_gen = S5PPRRR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 5, 6, 8, 9]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(7) = 0.0; % d5
