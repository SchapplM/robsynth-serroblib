% Umwandlung der Kinematikparameter von S4RPRR2V1 zu S4RPRR2
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4RPRR2V1
%   pkin_var=[a2 a3 d1 d3 theta2]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RPRR2
%   pkin_gen=[a2 a3 a4 d1 d3 d4 theta2]
%
% Siehe auch: S4RPRR2_structural_kinematic_parameters.m
function pkin_gen = S4RPRR2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([1, 2, 4, 5, 7]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(6) = 0.0; % d4
