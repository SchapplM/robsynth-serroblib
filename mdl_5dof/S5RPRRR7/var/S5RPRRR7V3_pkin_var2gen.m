% Umwandlung der Kinematikparameter von S5RPRRR7V3 zu S5RPRRR7
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RPRRR7V3
%   pkin_var=[a4 a5 d1 d4 d5 theta2]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RPRRR7
%   pkin_gen=[a2 a3 a4 a5 d1 d3 d4 d5 theta2]
%
% Siehe auch: S5RPRRR7_structural_kinematic_parameters.m
function pkin_gen = S5RPRRR7V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([3, 4, 5, 7, 8, 9]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(6) = 0.0; % d3
