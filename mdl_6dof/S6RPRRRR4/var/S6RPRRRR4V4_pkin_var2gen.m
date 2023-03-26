% Umwandlung der Kinematikparameter von S6RPRRRR4V4 zu S6RPRRRR4
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RPRRRR4V4
%   pkin_var=[a4 a5 d1 d4 d5 theta2]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RPRRRR4
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 d5 d6 theta2]
%
% Siehe auch: S6RPRRRR4_structural_kinematic_parameters.m
function pkin_gen = S6RPRRRR4V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([3, 4, 6, 8, 9, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % a6
pkin_gen(7) = 0.0; % d3
pkin_gen(10) = 0.0; % d6
