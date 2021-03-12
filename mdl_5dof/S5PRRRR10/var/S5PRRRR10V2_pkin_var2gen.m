% Umwandlung der Kinematikparameter von S5PRRRR10V2 zu S5PRRRR10
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5PRRRR10V2
%   pkin_var=[a2 a5 alpha2 d2 d5 theta1]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S5PRRRR10
%   pkin_gen=[a2 a3 a4 a5 alpha2 alpha3 d2 d3 d4 d5 theta1]
%
% Siehe auch: S5PRRRR10_structural_kinematic_parameters.m
function pkin_gen = S5PRRRR10V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 4, 5, 7, 10, 11]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(6) = pi/2; % alpha3
pkin_gen(8) = 0.0; % d3
pkin_gen(9) = 0.0; % d4
