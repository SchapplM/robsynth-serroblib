% Umwandlung der Kinematikparameter von S6RRRRPR4V3 zu S6RRRRPR4
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRRRPR4V3
%   pkin_var=[a2 a3 a5 a6 d1 d2 d3 d6 theta5]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRRPR4
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d6 theta5]
%
% Siehe auch: S6RRRRPR4_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR4V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 4, 5, 6, 7, 8, 10, 11]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(9) = 0.0; % d4
