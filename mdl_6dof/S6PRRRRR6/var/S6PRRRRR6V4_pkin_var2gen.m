% Umwandlung der Kinematikparameter von S6PRRRRR6V4 zu S6PRRRRR6
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6PRRRRR6V4
%   pkin_var=[a2 a4 alpha2 alpha4 d2 d4 theta1]
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6PRRRRR6
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d2 d3 d4 d5 d6 theta1]
%
% Siehe auch: S6PRRRRR6_structural_kinematic_parameters.m
function pkin_gen = S6PRRRRR6V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([1, 3, 6, 8, 9, 11, 14]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(7) = pi/2; % alpha3
pkin_gen(10) = 0.0; % d3
pkin_gen(12) = 0.0; % d5
pkin_gen(13) = 0.0; % d6
