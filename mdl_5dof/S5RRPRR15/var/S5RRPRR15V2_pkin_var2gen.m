% Umwandlung der Kinematikparameter von S5RRPRR15V2 zu S5RRPRR15
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRPRR15V2
%   pkin_var=[a5 d1 d5]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRPRR15
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 d5]
%
% Siehe auch: S5RRPRR15_structural_kinematic_parameters.m
function pkin_gen = S5RRPRR15V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([4, 5, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(6) = 0.0; % d2
pkin_gen(7) = 0.0; % d4
