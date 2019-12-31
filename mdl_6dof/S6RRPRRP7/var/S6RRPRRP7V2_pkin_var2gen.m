% Umwandlung der Kinematikparameter von S6RRPRRP7V2 zu S6RRPRRP7
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRPRRP7V2
%   pkin_var=[a2 a3 a4 a6 d1 d2 d4]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RRPRRP7
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4 d5]
%
% Siehe auch: S6RRPRRP7_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRP7V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 5, 6, 7, 8]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(9) = 0.0; % d5
