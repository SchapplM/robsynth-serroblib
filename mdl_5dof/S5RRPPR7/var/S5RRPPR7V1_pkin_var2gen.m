% Umwandlung der Kinematikparameter von S5RRPPR7V1 zu S5RRPPR7
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRPPR7V1
%   pkin_var=[a3 a4 a5 d1 d5 theta3]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRPPR7
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d5 theta3]
%
% Siehe auch: S5RRPPR7_structural_kinematic_parameters.m
function pkin_gen = S5RRPPR7V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([2, 3, 4, 5, 7, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(6) = 0.0; % d2
