% Umwandlung der Kinematikparameter von S6RRRRPR10V3 zu S6RRRRPR10
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRPR10V3
%   pkin_var=[a4 a5 a6 d1 d4 d6]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRRPR10
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d4 d6]
%
% Siehe auch: S6RRRRPR10_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR10V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([3, 4, 5, 7, 10, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(6) = pi/2; % alpha2
pkin_gen(8) = 0.0; % d2
pkin_gen(9) = 0.0; % d3
