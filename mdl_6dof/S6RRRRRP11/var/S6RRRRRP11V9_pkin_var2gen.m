% Umwandlung der Kinematikparameter von S6RRRRRP11V9 zu S6RRRRRP11
% Eingabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6RRRRRP11V9
%   pkin_var=[a2 a3 a4 a5 alpha2 alpha3 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRRP11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5]
%
% Siehe auch: S6RRRRRP11_structural_kinematic_parameters.m
function pkin_gen = S6RRRRRP11V9_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12]) = pkin_var;

pkin_gen(5) = 0.0; % a6
