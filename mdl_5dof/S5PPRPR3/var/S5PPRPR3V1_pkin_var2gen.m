% Umwandlung der Kinematikparameter von S5PPRPR3V1 zu S5PPRPR3
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S5PPRPR3V1
%   pkin_var=[a2 a3 a4 a5 d3 d5 theta1 theta4]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PPRPR3
%   pkin_gen=[a2 a3 a4 a5 d3 d5 theta1 theta2 theta4]
%
% Siehe auch: S5PPRPR3_structural_kinematic_parameters.m
function pkin_gen = S5PPRPR3V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 9]) = pkin_var;

pkin_gen(8) = pi/2; % theta2
