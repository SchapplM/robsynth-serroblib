% Umwandlung der Kinematikparameter von S6RRPRRR14V4 zu S6RRPRRR14
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14V4
%   pkin_var=[a3 a4 alpha3 alpha4 d1 d4 theta3]
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d1 d2 d4 d5 d6 theta3]
%
% Siehe auch: S6RRPRRR14_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRR14V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([2, 3, 7, 8, 9, 11, 14]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(10) = 0.0; % d2
pkin_gen(12) = 0.0; % d5
pkin_gen(13) = 0.0; % d6
