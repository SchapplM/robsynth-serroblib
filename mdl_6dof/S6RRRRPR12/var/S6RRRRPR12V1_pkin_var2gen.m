% Umwandlung der Kinematikparameter von S6RRRRPR12V1 zu S6RRRRPR12
% Eingabe:
% pkin_var (12x1) double
%   Kinematikparameter (pkin) von S6RRRRPR12V1
%   pkin_var=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 theta5]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRPR12
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d6 theta5]
%
% Siehe auch: S6RRRRPR12_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR12V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13]) = pkin_var;

pkin_gen(12) = 0.0; % d6
