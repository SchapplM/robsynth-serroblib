% Umwandlung der Kinematikparameter von S4RRRR6V1 zu S4RRRR6
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4RRRR6V1
%   pkin_var=[a2 alpha2 d1 d2]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S4RRRR6
%   pkin_gen=[a2 a3 a4 alpha2 d1 d2 d3 d4]
%
% Siehe auch: S4RRRR6_structural_kinematic_parameters.m
function pkin_gen = S4RRRR6V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 4, 5, 6]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(7) = 0.0; % d3
pkin_gen(8) = 0.0; % d4
