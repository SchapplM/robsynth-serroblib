% Umwandlung der Kinematikparameter von S4RRRR1V1 zu S4RRRR1
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4RRRR1V1
%   pkin_var=[a2 a3 d1 d2 d3]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RRRR1
%   pkin_gen=[a2 a3 a4 d1 d2 d3 d4]
%
% Siehe auch: S4RRRR1_structural_kinematic_parameters.m
function pkin_gen = S4RRRR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([1, 2, 4, 5, 6]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(7) = 0.0; % d4
