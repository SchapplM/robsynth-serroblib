% Umwandlung der Kinematikparameter von S4RPRR9V3 zu S4RPRR9
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RPRR9V3
%   pkin_var=[a4 d1 d4]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RPRR9
%   pkin_gen=[a2 a3 a4 d1 d3 d4]
%
% Siehe auch: S4RPRR9_structural_kinematic_parameters.m
function pkin_gen = S4RPRR9V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([3, 4, 6]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % d3
