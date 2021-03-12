% Umwandlung der Kinematikparameter von S6RRRRRR8V6 zu S6RRRRRR8
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8V6
%   pkin_var=[a4 a5 d1 d4 d5]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5 d6]
%
% Siehe auch: S6RRRRRR8_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRR8V6_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([3, 4, 8, 11, 12]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(7) = pi/2; % alpha3
pkin_gen(9) = 0.0; % d2
pkin_gen(10) = 0.0; % d3
pkin_gen(13) = 0.0; % d6
