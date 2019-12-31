% Umwandlung der Kinematikparameter von S6RRPRRR2V2 zu S6RRPRRR2
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRPRRR2V2
%   pkin_var=[a3 a4 a6 d1 d4 d6 theta3]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4 d5 d6 theta3]
%
% Siehe auch: S6RRPRRR2_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRR2V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([2, 3, 5, 6, 8, 10, 11]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(4) = 0.0; % a5
pkin_gen(7) = 0.0; % d2
pkin_gen(9) = 0.0; % d5
