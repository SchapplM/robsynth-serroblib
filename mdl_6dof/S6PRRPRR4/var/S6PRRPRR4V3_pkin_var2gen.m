% Umwandlung der Kinematikparameter von S6PRRPRR4V3 zu S6PRRPRR4
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRRPRR4V3
%   pkin_var=[a2 a3 a4 a5 alpha2 d2 d3 d5 theta1]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6PRRPRR4
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d3 d5 d6 theta1]
%
% Siehe auch: S6PRRPRR4_structural_kinematic_parameters.m
function pkin_gen = S6PRRPRR4V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 3, 4, 6, 7, 8, 9, 11]) = pkin_var;

pkin_gen(5) = 0.0; % a6
pkin_gen(10) = 0.0; % d6
