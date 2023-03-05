% Umwandlung der Kinematikparameter von S5RRRRP1V2 zu S5RRRRP1
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRRRP1V2
%   pkin_var=[a3 a4 d1 d3 d4]
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRRP1
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d4]
%
% Siehe auch: S5RRRRP1_structural_kinematic_parameters.m
function pkin_gen = S5RRRRP1V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([2, 3, 5, 7, 8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(4) = 0.0; % a5
pkin_gen(6) = 0.0; % d2
