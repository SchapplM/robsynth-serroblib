% Umwandlung der Kinematikparameter von S6RRRRRP11V2 zu S6RRRRRP11
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRRP11V2
%   pkin_var=[a2 a5 a6 alpha2 d1 d2 d5]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRRP11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5]
%
% Siehe auch: S6RRRRRP11_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRP11V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 4, 5, 6, 8, 9, 12]) = pkin_var;

pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(7) = pi/2; % alpha3
pkin_gen(10) = 0.0; % d3
pkin_gen(11) = 0.0; % d4
