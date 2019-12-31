% Umwandlung der Kinematikparameter von S4PRRR6V1 zu S4PRRR6
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4PRRR6V1
%   pkin_var=[a2 a4 d2 d4 theta1]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4PRRR6
%   pkin_gen=[a2 a3 a4 d2 d3 d4 theta1]
%
% Siehe auch: S4PRRR6_structural_kinematic_parameters.m
function pkin_gen = S4PRRR6V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([1, 3, 4, 6, 7]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % d3
