% Umwandlung der Kinematikparameter von S4RRRP2V3 zu S4RRRP2
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4RRRP2V3
%   pkin_var=[a2 a3 d1 d2 d3]
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RRRP2
%   pkin_gen=[a2 a3 a4 d1 d2 d3]
%
% Siehe auch: S4RRRP2_structural_kinematic_parameters.m
function pkin_gen = S4RRRP2V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([1, 2, 4, 5, 6]) = pkin_var;

pkin_gen(3) = 0.0; % a4
