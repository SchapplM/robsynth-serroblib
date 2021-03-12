% Umwandlung der Kinematikparameter von S4RRRR6V2 zu S4RRRR6
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRRR6V2
%   pkin_var=[a4 d1 d4]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S4RRRR6
%   pkin_gen=[a2 a3 a4 alpha2 d1 d2 d3 d4]
%
% Siehe auch: S4RRRR6_structural_kinematic_parameters.m
function pkin_gen = S4RRRR6V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([3, 5, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(4) = pi/2; % alpha2
pkin_gen(6) = 0.0; % d2
pkin_gen(7) = 0.0; % d3
