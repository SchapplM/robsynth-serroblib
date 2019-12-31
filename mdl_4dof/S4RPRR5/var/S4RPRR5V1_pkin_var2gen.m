% Umwandlung der Kinematikparameter von S4RPRR5V1 zu S4RPRR5
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4RPRR5V1
%   pkin_var=[a2 a3 d1 d3]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RPRR5
%   pkin_gen=[a2 a3 a4 d1 d3 d4]
%
% Siehe auch: S4RPRR5_structural_kinematic_parameters.m
function pkin_gen = S4RPRR5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([1, 2, 4, 5]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(6) = 0.0; % d4
