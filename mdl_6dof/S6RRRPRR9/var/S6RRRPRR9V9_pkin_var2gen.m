% Umwandlung der Kinematikparameter von S6RRRPRR9V9 zu S6RRRPRR9
% Eingabe:
% pkin_var (1x1) double
%   Kinematikparameter (pkin) von S6RRRPRR9V9
%   pkin_var=[d1]
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRPRR9
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d5 d6 theta4]
%
% Siehe auch: S6RRRPRR9_structural_kinematic_parameters.m
function pkin_gen = S6RRRPRR9V9_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([8]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(7) = pi/2; % alpha3
pkin_gen(9) = 0.0; % d2
pkin_gen(10) = 0.0; % d3
pkin_gen(11) = 0.0; % d5
pkin_gen(12) = 0.0; % d6
pkin_gen(13) = 0.0; % theta4
