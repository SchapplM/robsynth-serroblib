% Umwandlung der Kinematikparameter von S5PRRPR2V1 zu S5PRRPR2
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S5PRRPR2V1
%   pkin_var=[a2 a3 a4 a5 d2 d3 d5 theta1]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PRRPR2
%   pkin_gen=[a2 a3 a4 a5 d2 d3 d5 theta1 theta4]
%
% Siehe auch: S5PRRPR2_structural_kinematic_parameters.m
function pkin_gen = S5PRRPR2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8]) = pkin_var;

pkin_gen(9) = pi/2; % theta4
