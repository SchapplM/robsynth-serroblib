% Umwandlung der Kinematikparameter von S5RRPRR15V1 zu S5RRPRR15
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRPRR15V1
%   pkin_var=[a3 a4 a5 d1 d4 d5]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRPRR15
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 d5]
%
% Siehe auch: S5RRPRR15_structural_kinematic_parameters.m
function pkin_gen = S5RRPRR15V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([2, 3, 4, 5, 7, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(6) = 0.0; % d2
