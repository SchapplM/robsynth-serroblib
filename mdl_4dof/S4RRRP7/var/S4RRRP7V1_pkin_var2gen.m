% Umwandlung der Kinematikparameter von S4RRRP7V1 zu S4RRRP7
% Eingabe:
% pkin_var (2x1) double
%   Kinematikparameter (pkin) von S4RRRP7V1
%   pkin_var=[a4 d1]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RRRP7
%   pkin_gen=[a2 a3 a4 d1 d2 d3]
%
% Siehe auch: S4RRRP7_structural_kinematic_parameters.m
function pkin_gen = S4RRRP7V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([3, 4]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % d2
pkin_gen(6) = 0.0; % d3
