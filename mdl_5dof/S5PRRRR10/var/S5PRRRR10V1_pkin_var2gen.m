% Umwandlung der Kinematikparameter von S5PRRRR10V1 zu S5PRRRR10
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PRRRR10V1
%   pkin_var=[a2 a3 alpha2 alpha3 d2 d3 theta1]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S5PRRRR10
%   pkin_gen=[a2 a3 a4 a5 alpha2 alpha3 d2 d3 d4 d5 theta1]
%
% Siehe auch: S5PRRRR10_structural_kinematic_parameters.m
function pkin_gen = S5PRRRR10V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 5, 6, 7, 8, 11]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d4
pkin_gen(10) = 0.0; % d5
