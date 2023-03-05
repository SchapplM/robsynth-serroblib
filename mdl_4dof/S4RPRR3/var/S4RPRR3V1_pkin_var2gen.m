% Umwandlung der Kinematikparameter von S4RPRR3V1 zu S4RPRR3
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4RPRR3V1
%   pkin_var=[a4 d1 d4 theta2]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RPRR3
%   pkin_gen=[a2 a3 a4 d1 d3 d4 theta2]
%
% Siehe auch: S4RPRR3_structural_kinematic_parameters.m
function pkin_gen = S4RPRR3V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([3, 4, 6, 7]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % d3
