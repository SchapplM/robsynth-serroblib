% Umwandlung der Kinematikparameter von S5RPRRR12V1 zu S5RPRRR12
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RPRRR12V1
%   pkin_var=[a2 a3 a4 d1 d3 d4]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RPRRR12
%   pkin_gen=[a2 a3 a4 a5 d1 d3 d4 d5]
%
% Siehe auch: S5RPRRR12_structural_kinematic_parameters.m
function pkin_gen = S5RPRRR12V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 2, 3, 5, 6, 7]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(8) = 0.0; % d5
