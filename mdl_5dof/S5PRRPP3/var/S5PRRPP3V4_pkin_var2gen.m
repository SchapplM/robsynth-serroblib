% Umwandlung der Kinematikparameter von S5PRRPP3V4 zu S5PRRPP3
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5PRRPP3V4
%   pkin_var=[a2 a4 a5 d2 theta1]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5PRRPP3
%   pkin_gen=[a2 a3 a4 a5 d2 d3 theta1 theta4]
%
% Siehe auch: S5PRRPP3_structural_kinematic_parameters.m
function pkin_gen = S5PRRPP3V4_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 3, 4, 5, 7]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(6) = 0.0; % d3
pkin_gen(8) = 0.0; % theta4
