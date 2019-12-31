% Umwandlung der Kinematikparameter von S6PRPRRR4V2 zu S6PRPRRR4
% Eingabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6PRPRRR4V2
%   pkin_var=[a2 a3 a4 a6 alpha2 d2 d4 d6 theta1 theta3]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRPRRR4
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d4 d5 d6 theta1 theta3]
%
% Siehe auch: S6PRPRRR4_structural_kinematic_parameters.m
function pkin_gen = S6PRPRRR4V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 3, 5, 6, 7, 8, 10, 11, 12]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d5
