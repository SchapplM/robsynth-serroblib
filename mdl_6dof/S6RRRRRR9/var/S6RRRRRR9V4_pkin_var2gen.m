% Umwandlung der Kinematikparameter von S6RRRRRR9V4 zu S6RRRRRR9
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRRR9V4
%   pkin_var=[a3 a6 alpha3 d1 d3 d6]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRRR9
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5 d6]
%
% Siehe auch: S6RRRRRR9_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRR9V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([2, 5, 7, 8, 10, 13]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(6) = pi/2; % alpha2
pkin_gen(9) = 0.0; % d2
pkin_gen(11) = 0.0; % d4
pkin_gen(12) = 0.0; % d5
