% Umwandlung der Kinematikparameter von S6RRRPRR13V3 zu S6RRRPRR13
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRPRR13V3
%   pkin_var=[a2 a4 a5 alpha2 d1 d2 d5 theta4]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRPRR13
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d5 d6 theta4]
%
% Siehe auch: S6RRRPRR13_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRR13V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 3, 4, 6, 8, 9, 11, 13]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % a6
pkin_gen(7) = pi/2; % alpha3
pkin_gen(10) = 0.0; % d3
pkin_gen(12) = 0.0; % d6
