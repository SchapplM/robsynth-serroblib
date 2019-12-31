% Umwandlung der Kinematikparameter von S5RRPRR3V1 zu S5RRPRR3
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRPRR3V1
%   pkin_var=[a2 a3 a4 d1 d2 d4 theta3]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR3
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 d5 theta3]
%
% Siehe auch: S5RRPRR3_structural_kinematic_parameters.m
function pkin_gen = S5RRPRR3V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 5, 6, 7, 9]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(8) = 0.0; % d5
