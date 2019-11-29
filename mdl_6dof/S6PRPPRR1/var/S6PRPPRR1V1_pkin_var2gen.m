% Umwandlung der Kinematikparameter von S6PRPPRR1V1 zu S6PRPPRR1
% Eingabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6PRPPRR1V1
%   pkin_var=[a2 a3 a4 a5 a6 alpha2 d2 d5 theta1 theta3 theta4]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRPPRR1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d5 d6 theta1 theta3 theta4]
%
% Siehe auch: S6PRPPRR1_structural_kinematic_parameters.m
function pkin_gen = S6PRPPRR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12]) = pkin_var;

pkin_gen(9) = 0.0; % d6
