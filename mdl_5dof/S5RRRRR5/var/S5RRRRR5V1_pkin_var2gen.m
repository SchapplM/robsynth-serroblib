% Umwandlung der Kinematikparameter von S5RRRRR5V1 zu S5RRRRR5
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRRRR5V1
%   pkin_var=[a2 a3 a5 d1 d2 d3 d5]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRRR5
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4 d5]
%
% Siehe auch: S5RRRRR5_structural_kinematic_parameters.m
function pkin_gen = S5RRRRR5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 4, 5, 6, 7, 9]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(8) = 0.0; % d4
