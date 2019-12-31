% Umwandlung der Kinematikparameter von S5PRPRR9V1 zu S5PRPRR9
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5PRPRR9V1
%   pkin_var=[a2 a3 a4 d2 d4 theta1]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5PRPRR9
%   pkin_gen=[a2 a3 a4 a5 d2 d4 d5 theta1]
%
% Siehe auch: S5PRPRR9_structural_kinematic_parameters.m
function pkin_gen = S5PRPRR9V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 2, 3, 5, 6, 8]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(7) = 0.0; % d5
