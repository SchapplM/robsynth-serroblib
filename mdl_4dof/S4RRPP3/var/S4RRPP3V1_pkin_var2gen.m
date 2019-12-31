% Umwandlung der Kinematikparameter von S4RRPP3V1 zu S4RRPP3
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4RRPP3V1
%   pkin_var=[a3 a4 d1 theta3]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RRPP3
%   pkin_gen=[a2 a3 a4 d1 d2 theta3]
%
% Siehe auch: S4RRPP3_structural_kinematic_parameters.m
function pkin_gen = S4RRPP3V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([2, 3, 4, 6]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(5) = 0.0; % d2
