% Umwandlung der Kinematikparameter von S5RRRPP7V1 zu S5RRRPP7
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRRPP7V1
%   pkin_var=[a4 a5 d1]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S5RRRPP7
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3]
%
% Siehe auch: S5RRRPP7_structural_kinematic_parameters.m
function pkin_gen = S5RRRPP7V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([3, 4, 5]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(6) = 0.0; % d2
pkin_gen(7) = 0.0; % d3
