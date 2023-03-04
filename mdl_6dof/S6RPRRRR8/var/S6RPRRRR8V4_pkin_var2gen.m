% Umwandlung der Kinematikparameter von S6RPRRRR8V4 zu S6RPRRRR8
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RPRRRR8V4
%   pkin_var=[a4 a6 d1 d4 d6]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RPRRRR8
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 d5 d6]
%
% Siehe auch: S6RPRRRR8_structural_kinematic_parameters.m
function pkin_gen = S6RPRRRR8V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([3, 5, 6, 8, 10]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(4) = 0.0; % a5
pkin_gen(7) = 0.0; % d3
pkin_gen(9) = 0.0; % d5
