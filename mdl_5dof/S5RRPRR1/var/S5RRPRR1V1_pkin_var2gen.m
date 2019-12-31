% Umwandlung der Kinematikparameter von S5RRPRR1V1 zu S5RRPRR1
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRPRR1V1
%   pkin_var=[a3 a4 d4]
% Ausgabe:
% pkin_gen (4x1) double
%   Kinematikparameter (pkin) von S5RRPRR1
%   pkin_gen=[a3 a4 d4 d5]
%
% Siehe auch: S5RRPRR1_structural_kinematic_parameters.m
function pkin_gen = S5RRPRR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(4,1);
pkin_gen([1, 2, 3]) = pkin_var;

pkin_gen(4) = 0.0; % d5
