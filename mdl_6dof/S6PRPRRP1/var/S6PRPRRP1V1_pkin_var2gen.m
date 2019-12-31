% Umwandlung der Kinematikparameter von S6PRPRRP1V1 zu S6PRPRRP1
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRPRRP1V1
%   pkin_var=[a2 a3 a4 a6 alpha2 d2 d4 theta1 theta3]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6PRPRRP1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d4 d5 theta1 theta3]
%
% Siehe auch: S6PRPRRP1_structural_kinematic_parameters.m
function pkin_gen = S6PRPRRP1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 3, 5, 6, 7, 8, 10, 11]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d5
