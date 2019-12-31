% Umwandlung der Kinematikparameter von S5RRRRR1V1 zu S5RRRRR1
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRRRR1V1
%   pkin_var=[a2 a3 a4 d1]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S5RRRRR1
%   pkin_gen=[a2 a3 a4 a5 d1 d5]
%
% Siehe auch: S5RRRRR1_structural_kinematic_parameters.m
function pkin_gen = S5RRRRR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([1, 2, 3, 5]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(6) = 0.0; % d5
