% Umwandlung der Kinematikparameter von S6RRPRPR14V1 zu S6RRPRPR14
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRPRPR14V1
%   pkin_var=[a3 a4 a5 a6 d1 d4 d6]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPRPR14
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d6]
%
% Siehe auch: S6RRPRPR14_structural_kinematic_parameters.m
function pkin_gen = S6RRPRPR14V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([2, 3, 4, 5, 7, 9, 10]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(6) = pi/2; % alpha2
pkin_gen(8) = 0.0; % d2
