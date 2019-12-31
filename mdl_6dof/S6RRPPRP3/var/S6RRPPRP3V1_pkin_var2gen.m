% Umwandlung der Kinematikparameter von S6RRPPRP3V1 zu S6RRPPRP3
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPPRP3V1
%   pkin_var=[a3 a4 a5 a6 d1 d5]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S6RRPPRP3
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d5]
%
% Siehe auch: S6RRPPRP3_structural_kinematic_parameters.m
function pkin_gen = S6RRPPRP3V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([2, 3, 4, 5, 6, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(7) = 0.0; % d2
