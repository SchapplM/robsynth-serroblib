% Umwandlung der Kinematikparameter von S5RPRRR14V3 zu S5RPRRR14
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S5RPRRR14V3
%   pkin_var=[a4 a5 alpha2 alpha3 d1 d4 d5 theta2]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S5RPRRR14
%   pkin_gen=[a2 a3 a4 a5 alpha2 alpha3 d1 d3 d4 d5 theta2]
%
% Siehe auch: S5RPRRR14_structural_kinematic_parameters.m
function pkin_gen = S5RPRRR14V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([3, 4, 5, 6, 7, 9, 10, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(8) = 0.0; % d3
