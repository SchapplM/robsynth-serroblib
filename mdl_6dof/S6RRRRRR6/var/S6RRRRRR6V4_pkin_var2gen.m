% Umwandlung der Kinematikparameter von S6RRRRRR6V4 zu S6RRRRRR6
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RRRRRR6V4
%   pkin_var=[a4 a6 d1 d4 d6]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRRR6
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d4 d5 d6]
%
% Siehe auch: S6RRRRRR6_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRR6V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([3, 5, 7, 10, 12]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(4) = 0.0; % a5
pkin_gen(6) = pi/2; % alpha2
pkin_gen(8) = 0.0; % d2
pkin_gen(9) = 0.0; % d3
pkin_gen(11) = 0.0; % d5
