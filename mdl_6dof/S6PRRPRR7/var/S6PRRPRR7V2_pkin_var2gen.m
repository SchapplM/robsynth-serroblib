% Umwandlung der Kinematikparameter von S6PRRPRR7V2 zu S6PRRPRR7
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRRPRR7V2
%   pkin_var=[a2 a4 a5 a6 alpha2 d2 d5 d6 theta1]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6PRRPRR7
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d3 d5 d6 theta1]
%
% Siehe auch: S6PRRPRR7_structural_kinematic_parameters.m
function pkin_gen = S6PRRPRR7V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 3, 4, 5, 6, 7, 9, 10, 11]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(8) = 0.0; % d3
