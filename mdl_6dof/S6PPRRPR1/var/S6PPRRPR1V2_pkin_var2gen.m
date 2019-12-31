% Umwandlung der Kinematikparameter von S6PPRRPR1V2 zu S6PPRRPR1
% Eingabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6PPRRPR1V2
%   pkin_var=[a2 a3 a5 a6 alpha2 alpha3 d3 d6 theta1 theta2 theta5]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6PPRRPR1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d3 d4 d6 theta1 theta2 theta5]
%
% Siehe auch: S6PPRRPR1_structural_kinematic_parameters.m
function pkin_gen = S6PPRRPR1V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 2, 4, 5, 6, 7, 8, 10, 11, 12, 13]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(9) = 0.0; % d4
