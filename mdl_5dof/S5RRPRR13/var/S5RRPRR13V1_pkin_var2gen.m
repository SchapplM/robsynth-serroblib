% Umwandlung der Kinematikparameter von S5RRPRR13V1 zu S5RRPRR13
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRPRR13V1
%   pkin_var=[a3 a4 a5 d1 d4 d5 theta3]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR13
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 d5 theta3]
%
% Siehe auch: S5RRPRR13_structural_kinematic_parameters.m
function pkin_gen = S5RRPRR13V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([2, 3, 4, 5, 7, 8, 9]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(6) = 0.0; % d2
