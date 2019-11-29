% Umwandlung der Kinematikparameter von S6RPRRRR4V1 zu S6RPRRRR4
% Eingabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6RPRRRR4V1
%   pkin_var=[a2 a3 a4 a5 a6 d1 d3 d4 d5 theta2]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RPRRRR4
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 d5 d6 theta2]
%
% Siehe auch: S6RPRRRR4_structural_kinematic_parameters.m
function pkin_gen = S6RPRRRR4V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 9, 11]) = pkin_var;

pkin_gen(10) = 0.0; % d6
