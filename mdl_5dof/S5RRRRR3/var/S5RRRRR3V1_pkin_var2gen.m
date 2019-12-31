% Umwandlung der Kinematikparameter von S5RRRRR3V1 zu S5RRRRR3
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRRRR3V1
%   pkin_var=[a3 a5 d1]
% Ausgabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S5RRRRR3
%   pkin_gen=[a3 a4 a5 d1 d4]
%
% Siehe auch: S5RRRRR3_structural_kinematic_parameters.m
function pkin_gen = S5RRRRR3V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(5,1);
pkin_gen([1, 3, 4]) = pkin_var;

pkin_gen(2) = 0.0; % a4
pkin_gen(5) = 0.0; % d4
