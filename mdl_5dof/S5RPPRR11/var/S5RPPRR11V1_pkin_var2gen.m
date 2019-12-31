% Umwandlung der Kinematikparameter von S5RPPRR11V1 zu S5RPPRR11
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RPPRR11V1
%   pkin_var=[a2 a3 a4 d1 d4]
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S5RPPRR11
%   pkin_gen=[a2 a3 a4 a5 d1 d4 d5]
%
% Siehe auch: S5RPPRR11_structural_kinematic_parameters.m
function pkin_gen = S5RPPRR11V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([1, 2, 3, 5, 6]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(7) = 0.0; % d5
