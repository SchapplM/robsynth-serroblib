% Umwandlung der Kinematikparameter von S6RRPRRR12V3 zu S6RRPRRR12
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPRRR12V3
%   pkin_var=[a3 a4 a5 d1 d4 d5]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR12
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6]
%
% Siehe auch: S6RRPRRR12_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRR12V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([2, 3, 4, 7, 9, 10]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(8) = 0.0; % d2
pkin_gen(11) = 0.0; % d6
