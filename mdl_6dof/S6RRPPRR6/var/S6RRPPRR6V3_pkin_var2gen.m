% Umwandlung der Kinematikparameter von S6RRPPRR6V3 zu S6RRPPRR6
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRPPRR6V3
%   pkin_var=[a2 a3 a4 a5 d1 d2 d5 theta4]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPPRR6
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d5 d6 theta4]
%
% Siehe auch: S6RRPPRR6_structural_kinematic_parameters.m
function pkin_gen = S6RRPPRR6V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 4, 6, 7, 8, 10]) = pkin_var;

pkin_gen(5) = 0.0; % a6
pkin_gen(9) = 0.0; % d6
