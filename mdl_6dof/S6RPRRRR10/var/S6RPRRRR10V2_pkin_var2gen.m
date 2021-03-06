% Umwandlung der Kinematikparameter von S6RPRRRR10V2 zu S6RPRRRR10
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RPRRRR10V2
%   pkin_var=[a2 a3 a5 alpha2 alpha3 d1 d3 d5 theta2]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RPRRRR10
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d3 d4 d5 d6 theta2]
%
% Siehe auch: S6RPRRRR10_structural_kinematic_parameters.m
function pkin_gen = S6RPRRRR10V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 2, 4, 6, 7, 8, 9, 11, 13]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(5) = 0.0; % a6
pkin_gen(10) = 0.0; % d4
pkin_gen(12) = 0.0; % d6
