% Umwandlung der Kinematikparameter von S6RRRPRR8V3 zu S6RRRPRR8
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRPRR8V3
%   pkin_var=[a2 a6 alpha2 d1 d2 d6 theta4]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRPRR8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d5 d6 theta4]
%
% Siehe auch: S6RRRPRR8_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRR8V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 5, 6, 7, 8, 11, 12]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d3
pkin_gen(10) = 0.0; % d5
