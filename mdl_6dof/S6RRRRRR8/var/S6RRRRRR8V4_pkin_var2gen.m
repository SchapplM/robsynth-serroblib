% Umwandlung der Kinematikparameter von S6RRRRRR8V4 zu S6RRRRRR8
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8V4
%   pkin_var=[a2 a5 alpha2 d1 d2 d5]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5 d6]
%
% Siehe auch: S6RRRRRR8_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRR8V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 4, 6, 8, 9, 12]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(5) = 0.0; % a6
pkin_gen(7) = pi/2; % alpha3
pkin_gen(10) = 0.0; % d3
pkin_gen(11) = 0.0; % d4
pkin_gen(13) = 0.0; % d6
