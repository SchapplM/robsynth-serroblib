% Umwandlung der Kinematikparameter von S6RRRRRR10V4 zu S6RRRRRR10
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10V4
%   pkin_var=[a3 a4 d1 d4 d5 d6]
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d1 d2 d3 d4 d5 d6]
%
% Siehe auch: S6RRRRRR10_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRR10V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([2, 3, 9, 12, 13, 14]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(7) = 0.0; % alpha3
pkin_gen(8) = 0.0; % alpha4
pkin_gen(10) = 0.0; % d2
pkin_gen(11) = 0.0; % d3
