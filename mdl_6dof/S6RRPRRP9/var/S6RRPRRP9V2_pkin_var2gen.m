% Umwandlung der Kinematikparameter von S6RRPRRP9V2 zu S6RRPRRP9
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPRRP9V2
%   pkin_var=[a3 a4 a6 d1 d4 theta3]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRP9
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 theta3]
%
% Siehe auch: S6RRPRRP9_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRP9V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([2, 3, 5, 7, 9, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(4) = 0.0; % a5
pkin_gen(6) = pi/2; % alpha2
pkin_gen(8) = 0.0; % d2
pkin_gen(10) = 0.0; % d5
