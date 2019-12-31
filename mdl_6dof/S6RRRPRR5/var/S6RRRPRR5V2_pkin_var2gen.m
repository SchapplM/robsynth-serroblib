% Umwandlung der Kinematikparameter von S6RRRPRR5V2 zu S6RRRPRR5
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRPRR5V2
%   pkin_var=[a3 a4 a5 a6 d1 d3 d5 d6]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRPRR5
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d5 d6]
%
% Siehe auch: S6RRRPRR5_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRR5V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([2, 3, 4, 5, 6, 8, 9, 10]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(7) = 0.0; % d2
