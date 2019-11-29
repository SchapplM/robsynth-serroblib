% Umwandlung der Kinematikparameter von S6RPRRPR12V1 zu S6RPRRPR12
% Eingabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6RPRRPR12V1
%   pkin_var=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d3 d4 theta2]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RPRRPR12
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d3 d4 d6 theta2]
%
% Siehe auch: S6RPRRPR12_structural_kinematic_parameters.m
function pkin_gen = S6RPRRPR12V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]) = pkin_var;

pkin_gen(11) = 0.0; % d6
