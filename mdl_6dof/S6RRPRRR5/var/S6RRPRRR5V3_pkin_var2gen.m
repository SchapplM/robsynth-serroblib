% Umwandlung der Kinematikparameter von S6RRPRRR5V3 zu S6RRPRRR5
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRPRRR5V3
%   pkin_var=[a3 a4 a6 d1 d4 d6 theta3]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRPRRR5
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6 theta3]
%
% Siehe auch: S6RRPRRR5_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRR5V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([2, 3, 5, 7, 9, 11, 12]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(4) = 0.0; % a5
pkin_gen(6) = pi/2; % alpha2
pkin_gen(8) = 0.0; % d2
pkin_gen(10) = 0.0; % d5
