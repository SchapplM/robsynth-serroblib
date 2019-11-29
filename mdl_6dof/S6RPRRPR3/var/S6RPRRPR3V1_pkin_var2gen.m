% Umwandlung der Kinematikparameter von S6RPRRPR3V1 zu S6RPRRPR3
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RPRRPR3V1
%   pkin_var=[a2 a3 a4 a5 a6 d1 d3 d4 theta2]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RPRRPR3
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 d6 theta2]
%
% Siehe auch: S6RPRRPR3_structural_kinematic_parameters.m
function pkin_gen = S6RPRRPR3V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 10]) = pkin_var;

pkin_gen(9) = 0.0; % d6
