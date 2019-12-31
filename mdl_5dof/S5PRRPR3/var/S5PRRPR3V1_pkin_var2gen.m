% Umwandlung der Kinematikparameter von S5PRRPR3V1 zu S5PRRPR3
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PRRPR3V1
%   pkin_var=[a2 a4 a5 d2 d5 theta1 theta4]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PRRPR3
%   pkin_gen=[a2 a3 a4 a5 d2 d3 d5 theta1 theta4]
%
% Siehe auch: S5PRRPR3_structural_kinematic_parameters.m
function pkin_gen = S5PRRPR3V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 3, 4, 5, 7, 8, 9]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(6) = 0.0; % d3
