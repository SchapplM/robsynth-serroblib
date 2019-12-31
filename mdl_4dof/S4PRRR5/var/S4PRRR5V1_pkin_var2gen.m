% Umwandlung der Kinematikparameter von S4PRRR5V1 zu S4PRRR5
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4PRRR5V1
%   pkin_var=[a2 a3 d2 d3 theta1]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4PRRR5
%   pkin_gen=[a2 a3 a4 d2 d3 d4 theta1]
%
% Siehe auch: S4PRRR5_structural_kinematic_parameters.m
function pkin_gen = S4PRRR5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([1, 2, 4, 5, 7]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(6) = 0.0; % d4
