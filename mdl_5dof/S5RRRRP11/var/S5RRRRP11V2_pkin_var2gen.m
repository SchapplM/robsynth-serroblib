% Umwandlung der Kinematikparameter von S5RRRRP11V2 zu S5RRRRP11
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRRRP11V2
%   pkin_var=[a4 a5 d1 d4]
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRRP11
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d4]
%
% Siehe auch: S5RRRRP11_structural_kinematic_parameters.m
function pkin_gen = S5RRRRP11V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([3, 4, 6, 9]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(5) = pi/2; % alpha2
pkin_gen(7) = 0.0; % d2
pkin_gen(8) = 0.0; % d3
