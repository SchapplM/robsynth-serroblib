% Umwandlung der Kinematikparameter von S6RRRPRR1V2 zu S6RRRPRR1
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRPRR1V2
%   pkin_var=[a3 a4 a5 d1 d3 d5 theta4]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRPRR1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d5 d6 theta4]
%
% Siehe auch: S6RRRPRR1_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRR1V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([2, 3, 4, 6, 8, 9, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(5) = 0.0; % a6
pkin_gen(7) = 0.0; % d2
pkin_gen(10) = 0.0; % d6
