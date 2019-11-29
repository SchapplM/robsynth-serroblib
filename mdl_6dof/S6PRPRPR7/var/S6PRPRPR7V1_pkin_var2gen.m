% Umwandlung der Kinematikparameter von S6PRPRPR7V1 zu S6PRPRPR7
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRPRPR7V1
%   pkin_var=[a2 a3 a4 a5 a6 alpha2 d2 d4 theta1]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6PRPRPR7
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d4 d6 theta1]
%
% Siehe auch: S6PRPRPR7_structural_kinematic_parameters.m
function pkin_gen = S6PRPRPR7V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 10]) = pkin_var;

pkin_gen(9) = 0.0; % d6
