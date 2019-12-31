% Umwandlung der Kinematikparameter von S6PPRRRP2V1 zu S6PPRRRP2
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6PPRRRP2V1
%   pkin_var=[a2 a3 a6 alpha2 alpha3 d3 theta1 theta2]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PPRRRP2
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d3 d4 d5 theta1 theta2]
%
% Siehe auch: S6PPRRRP2_structural_kinematic_parameters.m
function pkin_gen = S6PPRRRP2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 5, 6, 7, 8, 11, 12]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d4
pkin_gen(10) = 0.0; % d5
