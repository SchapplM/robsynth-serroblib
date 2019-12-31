% Umwandlung der Kinematikparameter von S6RRPRRP2V1 zu S6RRPRRP2
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPRRP2V1
%   pkin_var=[a3 a4 a6 d1 d4 theta3]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPRRP2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4 d5 theta3]
%
% Siehe auch: S6RRPRRP2_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRP2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([2, 3, 5, 6, 8, 10]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(4) = 0.0; % a5
pkin_gen(7) = 0.0; % d2
pkin_gen(9) = 0.0; % d5
