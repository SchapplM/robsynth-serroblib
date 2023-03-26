% Umwandlung der Kinematikparameter von S6RRRPRR9V10 zu S6RRRPRR9
% Eingabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6RRRPRR9V10
%   pkin_var=[a2 a3 a6 alpha2 alpha3 d1 d2 d3 d6 theta4]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRPRR9
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d5 d6 theta4]
%
% Siehe auch: S6RRRPRR9_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRR9V10_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 2, 5, 6, 7, 8, 9, 10, 12, 13]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(11) = 0.0; % d5
