% Umwandlung der Kinematikparameter von S5RRPPR6V1 zu S5RRPPR6
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRPPR6V1
%   pkin_var=[a3 a4 a5 d1 d5 theta3 theta4]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPPR6
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d5 theta3 theta4]
%
% Siehe auch: S5RRPPR6_structural_kinematic_parameters.m
function pkin_gen = S5RRPPR6V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([2, 3, 4, 5, 7, 8, 9]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(6) = 0.0; % d2
