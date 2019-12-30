% Umwandlung der Kinematikparameter von S5PRPRP2V1 zu S5PRPRP2
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PRPRP2V1
%   pkin_var=[a2 a3 a4 a5 d2 d4 theta1]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5PRPRP2
%   pkin_gen=[a2 a3 a4 a5 d2 d4 theta1 theta3]
%
% Siehe auch: S5PRPRP2_structural_kinematic_parameters.m
function pkin_gen = S5PRPRP2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7]) = pkin_var;

pkin_gen(8) = pi/2; % theta3
