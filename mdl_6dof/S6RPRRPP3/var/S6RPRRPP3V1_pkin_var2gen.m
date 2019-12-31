% Umwandlung der Kinematikparameter von S6RPRRPP3V1 zu S6RPRRPP3
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPRRPP3V1
%   pkin_var=[a2 a3 a5 a6 d1 d3 theta2]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPRRPP3
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 theta2]
%
% Siehe auch: S6RPRRPP3_structural_kinematic_parameters.m
function pkin_gen = S6RPRRPP3V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 4, 5, 6, 7, 9]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(8) = 0.0; % d4
