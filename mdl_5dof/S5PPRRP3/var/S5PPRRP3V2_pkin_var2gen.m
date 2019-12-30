% Umwandlung der Kinematikparameter von S5PPRRP3V2 zu S5PPRRP3
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PPRRP3V2
%   pkin_var=[a2 a3 a4 a5 d3 d4 theta1]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5PPRRP3
%   pkin_gen=[a2 a3 a4 a5 d3 d4 theta1 theta2]
%
% Siehe auch: S5PPRRP3_structural_kinematic_parameters.m
function pkin_gen = S5PPRRP3V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7]) = pkin_var;

pkin_gen(8) = pi/2; % theta2
