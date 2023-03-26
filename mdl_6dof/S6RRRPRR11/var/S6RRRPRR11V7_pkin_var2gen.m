% Umwandlung der Kinematikparameter von S6RRRPRR11V7 zu S6RRRPRR11
% Eingabe:
% pkin_var (1x1) double
%   Kinematikparameter (pkin) von S6RRRPRR11V7
%   pkin_var=[d1]
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRPRR11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d5 d6]
%
% Siehe auch: S6RRRPRR11_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRR11V7_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([7]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(8) = 0.0; % d2
pkin_gen(9) = 0.0; % d3
pkin_gen(10) = 0.0; % d5
pkin_gen(11) = 0.0; % d6
