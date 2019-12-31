% Umwandlung der Kinematikparameter von S4RRRR5V1 zu S4RRRR5
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRRR5V1
%   pkin_var=[a4 d1 d4]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RRRR5
%   pkin_gen=[a2 a3 a4 d1 d2 d3 d4]
%
% Siehe auch: S4RRRR5_structural_kinematic_parameters.m
function pkin_gen = S4RRRR5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([3, 4, 7]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(5) = 0.0; % d2
pkin_gen(6) = 0.0; % d3
