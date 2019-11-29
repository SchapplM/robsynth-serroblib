% Umwandlung der Kinematikparameter von S4PPPR1V1 zu S4PPPR1
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4PPPR1V1
%   pkin_var=[a2 a3 a4 theta1]
% Ausgabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S4PPPR1
%   pkin_gen=[a2 a3 a4 d4 theta1]
%
% Siehe auch: S4PPPR1_structural_kinematic_parameters.m
function pkin_gen = S4PPPR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(5,1);
pkin_gen([1, 2, 3, 5]) = pkin_var;

pkin_gen(4) = 0.0; % d4
