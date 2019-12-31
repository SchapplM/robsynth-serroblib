% Umwandlung der Kinematikparameter von S6RPPPRR1V2 zu S6RPPPRR1
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPPPRR1V2
%   pkin_var=[a2 a3 a4 a5 d1 d5 theta2]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPPPRR1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d5 d6 theta2]
%
% Siehe auch: S6RPPPRR1_structural_kinematic_parameters.m
function pkin_gen = S6RPPPRR1V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 4, 6, 7, 9]) = pkin_var;

pkin_gen(5) = 0.0; % a6
pkin_gen(8) = 0.0; % d6
