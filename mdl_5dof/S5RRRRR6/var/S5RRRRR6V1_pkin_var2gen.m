% Umwandlung der Kinematikparameter von S5RRRRR6V1 zu S5RRRRR6
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRRRR6V1
%   pkin_var=[a2 a4 a5 d1 d2 d4 d5]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRRR6
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4 d5]
%
% Siehe auch: S5RRRRR6_structural_kinematic_parameters.m
function pkin_gen = S5RRRRR6V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 3, 4, 5, 6, 8, 9]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(7) = 0.0; % d3
