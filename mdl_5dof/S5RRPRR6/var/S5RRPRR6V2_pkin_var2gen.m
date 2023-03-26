% Umwandlung der Kinematikparameter von S5RRPRR6V2 zu S5RRPRR6
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRPRR6V2
%   pkin_var=[a2 a5 d1 d2 d5]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR6
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 d5 theta3]
%
% Siehe auch: S5RRPRR6_structural_kinematic_parameters.m
function pkin_gen = S5RRPRR6V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 4, 5, 6, 8]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(7) = 0.0; % d4
pkin_gen(9) = pi/2; % theta3
