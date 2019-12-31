% Umwandlung der Kinematikparameter von S6RRPRPP5V1 zu S6RRPRPP5
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPRPP5V1
%   pkin_var=[a3 a4 a5 a6 d1 d4]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S6RRPRPP5
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4]
%
% Siehe auch: S6RRPRPP5_structural_kinematic_parameters.m
function pkin_gen = S6RRPRPP5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([2, 3, 4, 5, 6, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(7) = 0.0; % d2
