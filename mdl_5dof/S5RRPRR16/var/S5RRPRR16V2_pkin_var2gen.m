% Umwandlung der Kinematikparameter von S5RRPRR16V2 zu S5RRPRR16
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRPRR16V2
%   pkin_var=[a3 a4 d1 d4]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR16
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d4 d5]
%
% Siehe auch: S5RRPRR16_structural_kinematic_parameters.m
function pkin_gen = S5RRPRR16V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([2, 3, 6, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = pi/2; % alpha2
pkin_gen(7) = 0.0; % d2
pkin_gen(9) = 0.0; % d5
