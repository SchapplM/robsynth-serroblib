% Umwandlung der Kinematikparameter von S4PRRP5V1 zu S4PRRP5
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4PRRP5V1
%   pkin_var=[a2 a4 d2 theta1]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4PRRP5
%   pkin_gen=[a2 a3 a4 d2 d3 theta1]
%
% Siehe auch: S4PRRP5_structural_kinematic_parameters.m
function pkin_gen = S4PRRP5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([1, 3, 4, 6]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % d3
