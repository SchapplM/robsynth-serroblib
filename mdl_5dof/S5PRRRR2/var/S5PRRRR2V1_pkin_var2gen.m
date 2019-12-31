% Umwandlung der Kinematikparameter von S5PRRRR2V1 zu S5PRRRR2
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5PRRRR2V1
%   pkin_var=[a2 a3 a4 d3 d4]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S5PRRRR2
%   pkin_gen=[a2 a3 a4 d3 d4 d5]
%
% Siehe auch: S5PRRRR2_structural_kinematic_parameters.m
function pkin_gen = S5PRRRR2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([1, 2, 3, 4, 5]) = pkin_var;

pkin_gen(6) = 0.0; % d5
