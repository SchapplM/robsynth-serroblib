% Umwandlung der Kinematikparameter von S5RRPRR14V5 zu S5RRPRR14
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRPRR14V5
%   pkin_var=[a2 a5 alpha2 d1 d2 d5 theta3]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5RRPRR14
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d4 d5 theta3]
%
% Siehe auch: S5RRPRR14_structural_kinematic_parameters.m
function pkin_gen = S5RRPRR14V5_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 4, 5, 6, 7, 9, 10]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(8) = 0.0; % d4
