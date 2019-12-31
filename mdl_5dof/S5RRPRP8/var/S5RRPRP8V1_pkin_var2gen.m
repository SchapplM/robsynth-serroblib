% Umwandlung der Kinematikparameter von S5RRPRP8V1 zu S5RRPRP8
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRPRP8V1
%   pkin_var=[a3 a4 a5 d1 d4]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S5RRPRP8
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4]
%
% Siehe auch: S5RRPRP8_structural_kinematic_parameters.m
function pkin_gen = S5RRPRP8V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([2, 3, 4, 5, 7]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(6) = 0.0; % d2
