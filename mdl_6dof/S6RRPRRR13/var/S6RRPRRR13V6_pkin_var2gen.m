% Umwandlung der Kinematikparameter von S6RRPRRR13V6 zu S6RRPRRR13
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13V6
%   pkin_var=[a2 a5 a6 alpha2 d1 d2 d5 d6]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR13
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d5 d6]
%
% Siehe auch: S6RRPRRR13_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRR13V6_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 4, 5, 6, 7, 8, 10, 11]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(9) = 0.0; % d4
