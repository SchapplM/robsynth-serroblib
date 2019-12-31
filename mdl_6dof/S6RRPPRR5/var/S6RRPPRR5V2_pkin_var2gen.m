% Umwandlung der Kinematikparameter von S6RRPPRR5V2 zu S6RRPPRR5
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRPPRR5V2
%   pkin_var=[a2 a3 a4 a5 alpha2 d1 d2 d5]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPPRR5
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d5 d6]
%
% Siehe auch: S6RRPPRR5_structural_kinematic_parameters.m
function pkin_gen = S6RRPPRR5V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 4, 6, 7, 8, 9]) = pkin_var;

pkin_gen(5) = 0.0; % a6
pkin_gen(10) = 0.0; % d6
