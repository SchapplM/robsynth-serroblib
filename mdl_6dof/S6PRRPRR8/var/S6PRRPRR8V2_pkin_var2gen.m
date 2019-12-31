% Umwandlung der Kinematikparameter von S6PRRPRR8V2 zu S6PRRPRR8
% Eingabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6PRRPRR8V2
%   pkin_var=[a2 a3 a4 a5 alpha2 alpha3 d2 d3 d5 theta1]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRRPRR8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d2 d3 d5 d6 theta1]
%
% Siehe auch: S6PRRPRR8_structural_kinematic_parameters.m
function pkin_gen = S6PRRPRR8V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 3, 4, 6, 7, 8, 9, 10, 12]) = pkin_var;

pkin_gen(5) = 0.0; % a6
pkin_gen(11) = 0.0; % d6
