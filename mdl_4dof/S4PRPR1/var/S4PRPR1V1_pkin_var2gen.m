% Umwandlung der Kinematikparameter von S4PRPR1V1 zu S4PRPR1
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4PRPR1V1
%   pkin_var=[a2 a3 a4 d2 theta1]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4PRPR1
%   pkin_gen=[a2 a3 a4 d2 d4 theta1]
%
% Siehe auch: S4PRPR1_structural_kinematic_parameters.m
function pkin_gen = S4PRPR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([1, 2, 3, 4, 6]) = pkin_var;

pkin_gen(5) = 0.0; % d4
