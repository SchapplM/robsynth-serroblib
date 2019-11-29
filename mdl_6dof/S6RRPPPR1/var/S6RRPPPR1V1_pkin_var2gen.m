% Umwandlung der Kinematikparameter von S6RRPPPR1V1 zu S6RRPPPR1
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRPPPR1V1
%   pkin_var=[a2 a3 a4 a5 a6 d1 d2 theta3 theta4]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPPPR1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d6 theta3 theta4]
%
% Siehe auch: S6RRPPPR1_structural_kinematic_parameters.m
function pkin_gen = S6RRPPPR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 9, 10]) = pkin_var;

pkin_gen(8) = 0.0; % d6
