% Umwandlung der Kinematikparameter von S6PRRPRR8V3 zu S6PRRPRR8
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6PRRPRR8V3
%   pkin_var=[a2 a4 a5 alpha2 d2 d5 theta1]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRRPRR8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d2 d3 d5 d6 theta1]
%
% Siehe auch: S6PRRPRR8_structural_kinematic_parameters.m
function pkin_gen = S6PRRPRR8V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 3, 4, 6, 8, 10, 12]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % a6
pkin_gen(7) = pi/2; % alpha3
pkin_gen(9) = 0.0; % d3
pkin_gen(11) = 0.0; % d6
