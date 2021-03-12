% Umwandlung der Kinematikparameter von S6PRRRRR6V3 zu S6PRRRRR6
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRRRRR6V3
%   pkin_var=[a2 a3 a6 alpha2 alpha3 d2 d3 d6 theta1]
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6PRRRRR6
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d2 d3 d4 d5 d6 theta1]
%
% Siehe auch: S6PRRRRR6_structural_kinematic_parameters.m
function pkin_gen = S6PRRRRR6V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([1, 2, 5, 6, 7, 9, 10, 13, 14]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(8) = pi/2; % alpha4
pkin_gen(11) = 0.0; % d4
pkin_gen(12) = 0.0; % d5
