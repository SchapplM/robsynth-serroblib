% Umwandlung der Kinematikparameter von S4RRPR8V3 zu S4RRPR8
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRPR8V3
%   pkin_var=[a2 d1 d2]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RRPR8
%   pkin_gen=[a2 a3 a4 d1 d2 d4]
%
% Siehe auch: S4RRPR8_structural_kinematic_parameters.m
function pkin_gen = S4RRPR8V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([1, 4, 5]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(6) = 0.0; % d4
