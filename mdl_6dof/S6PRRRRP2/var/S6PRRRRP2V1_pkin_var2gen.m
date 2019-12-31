% Umwandlung der Kinematikparameter von S6PRRRRP2V1 zu S6PRRRRP2
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6PRRRRP2V1
%   pkin_var=[a2 a4 a6 alpha2 d2 d4 theta1]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6PRRRRP2
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d3 d4 d5 theta1]
%
% Siehe auch: S6PRRRRP2_structural_kinematic_parameters.m
function pkin_gen = S6PRRRRP2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 3, 5, 6, 7, 9, 11]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(4) = 0.0; % a5
pkin_gen(8) = 0.0; % d3
pkin_gen(10) = 0.0; % d5
