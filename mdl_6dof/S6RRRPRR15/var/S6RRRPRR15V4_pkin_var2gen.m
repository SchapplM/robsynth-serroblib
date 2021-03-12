% Umwandlung der Kinematikparameter von S6RRRPRR15V4 zu S6RRRPRR15
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRPRR15V4
%   pkin_var=[a3 a4 a5 alpha3 d1 d3 d5]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRPRR15
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d5 d6]
%
% Siehe auch: S6RRRPRR15_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRR15V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([2, 3, 4, 7, 8, 10, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(9) = 0.0; % d2
pkin_gen(12) = 0.0; % d6
