% Umwandlung der Kinematikparameter von S5PRRRP5V1 zu S5PRRRP5
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5PRRRP5V1
%   pkin_var=[a2 a4 a5 d2 d4 theta1]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5PRRRP5
%   pkin_gen=[a2 a3 a4 a5 d2 d3 d4 theta1]
%
% Siehe auch: S5PRRRP5_structural_kinematic_parameters.m
function pkin_gen = S5PRRRP5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 3, 4, 5, 7, 8]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(6) = 0.0; % d3
