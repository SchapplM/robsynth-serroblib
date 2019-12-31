% Umwandlung der Kinematikparameter von S5PRRRR8V2 zu S5PRRRR8
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S5PRRRR8V2
%   pkin_var=[a2 a3 a4 alpha2 d2 d3 d4 theta1]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5PRRRR8
%   pkin_gen=[a2 a3 a4 a5 alpha2 d2 d3 d4 d5 theta1]
%
% Siehe auch: S5PRRRR8_structural_kinematic_parameters.m
function pkin_gen = S5PRRRR8V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 5, 6, 7, 8, 10]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d5
