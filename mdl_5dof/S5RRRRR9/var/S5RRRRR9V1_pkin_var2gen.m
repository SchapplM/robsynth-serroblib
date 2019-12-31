% Umwandlung der Kinematikparameter von S5RRRRR9V1 zu S5RRRRR9
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRRRR9V1
%   pkin_var=[a4 a5 d1 d4 d5]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRRR9
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4 d5]
%
% Siehe auch: S5RRRRR9_structural_kinematic_parameters.m
function pkin_gen = S5RRRRR9V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([3, 4, 5, 8, 9]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(6) = 0.0; % d2
pkin_gen(7) = 0.0; % d3
