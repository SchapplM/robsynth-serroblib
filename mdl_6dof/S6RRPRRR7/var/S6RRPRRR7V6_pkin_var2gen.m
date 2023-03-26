% Umwandlung der Kinematikparameter von S6RRPRRR7V6 zu S6RRPRRR7
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RRPRRR7V6
%   pkin_var=[a2 a6 d1 d2 d6]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPRRR7
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4 d5 d6]
%
% Siehe auch: S6RRPRRR7_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRR7V6_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 5, 6, 7, 10]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(8) = 0.0; % d4
pkin_gen(9) = 0.0; % d5
