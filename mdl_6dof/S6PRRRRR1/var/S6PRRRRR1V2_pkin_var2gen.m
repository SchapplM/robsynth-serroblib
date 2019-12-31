% Umwandlung der Kinematikparameter von S6PRRRRR1V2 zu S6PRRRRR1
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6PRRRRR1V2
%   pkin_var=[a2 a4 a5 alpha2 d2 d4 d5 theta1]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRRRRR1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d3 d4 d5 d6 theta1]
%
% Siehe auch: S6PRRRRR1_structural_kinematic_parameters.m
function pkin_gen = S6PRRRRR1V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 3, 4, 6, 7, 9, 10, 12]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % a6
pkin_gen(8) = 0.0; % d3
pkin_gen(11) = 0.0; % d6
