% Umwandlung der Kinematikparameter von S5RPRRP2V1 zu S5RPRRP2
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RPRRP2V1
%   pkin_var=[a2 a3 a5 d1 d3 theta2]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RPRRP2
%   pkin_gen=[a2 a3 a4 a5 d1 d3 d4 theta2]
%
% Siehe auch: S5RPRRP2_structural_kinematic_parameters.m
function pkin_gen = S5RPRRP2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 2, 4, 5, 6, 8]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(7) = 0.0; % d4
