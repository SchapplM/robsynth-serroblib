% Umwandlung der Kinematikparameter von S4RRPP5V1 zu S4RRPP5
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRPP5V1
%   pkin_var=[a3 a4 d1]
% Ausgabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S4RRPP5
%   pkin_gen=[a2 a3 a4 d1 d2]
%
% Siehe auch: S4RRPP5_structural_kinematic_parameters.m
function pkin_gen = S4RRPP5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(5,1);
pkin_gen([2, 3, 4]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(5) = 0.0; % d2
