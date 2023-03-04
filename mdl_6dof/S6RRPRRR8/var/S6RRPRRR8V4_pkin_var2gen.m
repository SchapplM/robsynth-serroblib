% Umwandlung der Kinematikparameter von S6RRPRRR8V4 zu S6RRPRRR8
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPRRR8V4
%   pkin_var=[a5 a6 d1 d5 d6 theta3]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR8
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4 d5 d6 theta3]
%
% Siehe auch: S6RRPRRR8_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRR8V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([4, 5, 6, 9, 10, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(7) = 0.0; % d2
pkin_gen(8) = 0.0; % d4
