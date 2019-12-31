% Umwandlung der Kinematikparameter von S4PPRR4V1 zu S4PPRR4
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4PPRR4V1
%   pkin_var=[a2 a3 d3 theta1 theta2]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4PPRR4
%   pkin_gen=[a2 a3 a4 d3 d4 theta1 theta2]
%
% Siehe auch: S4PPRR4_structural_kinematic_parameters.m
function pkin_gen = S4PPRR4V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([1, 2, 4, 6, 7]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(5) = 0.0; % d4
