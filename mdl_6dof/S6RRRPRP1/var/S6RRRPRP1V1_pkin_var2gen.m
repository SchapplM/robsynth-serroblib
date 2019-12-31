% Umwandlung der Kinematikparameter von S6RRRPRP1V1 zu S6RRRPRP1
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRPRP1V1
%   pkin_var=[a3 a4 a5 a6 d1 d3 d5 theta4]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRPRP1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d5 theta4]
%
% Siehe auch: S6RRRPRP1_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRP1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([2, 3, 4, 5, 6, 8, 9, 10]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(7) = 0.0; % d2
