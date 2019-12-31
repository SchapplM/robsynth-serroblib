% Umwandlung der Kinematikparameter von S5RRRPP5V1 zu S5RRRPP5
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRRPP5V1
%   pkin_var=[a3 a4 a5 d1 d3]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S5RRRPP5
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3]
%
% Siehe auch: S5RRRPP5_structural_kinematic_parameters.m
function pkin_gen = S5RRRPP5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([2, 3, 4, 5, 7]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(6) = 0.0; % d2
