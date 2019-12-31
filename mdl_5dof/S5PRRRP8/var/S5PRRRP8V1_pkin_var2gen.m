% Umwandlung der Kinematikparameter von S5PRRRP8V1 zu S5PRRRP8
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5PRRRP8V1
%   pkin_var=[a2 a5 alpha2 d2 theta1]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PRRRP8
%   pkin_gen=[a2 a3 a4 a5 alpha2 d2 d3 d4 theta1]
%
% Siehe auch: S5PRRRP8_structural_kinematic_parameters.m
function pkin_gen = S5PRRRP8V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 4, 5, 6, 9]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(7) = 0.0; % d3
pkin_gen(8) = 0.0; % d4
