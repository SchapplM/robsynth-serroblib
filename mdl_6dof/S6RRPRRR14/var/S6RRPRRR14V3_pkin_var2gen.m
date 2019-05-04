% Umwandlung der Kinematikparameter von S6RRPRRR14V3 zu S6RRPRRR14
% Eingabe:
% pkin_var (0x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14V3
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14
%
% Siehe auch: S6RRPRRR14_structural_kinematic_parameters.m
function pkin_gen = S6RRPRRR14V3_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(2) = 0.0; % a3
pkin_gen(3) = 0.0; % a4
pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(6) = pi/2; % alpha2
pkin_gen(7) = pi/2; % alpha3
pkin_gen(8) = 0.0; % alpha4
pkin_gen(9) = 0.0; % d1
pkin_gen(10) = 0.0; % d2
pkin_gen(11) = 0.0; % d4
pkin_gen(12) = 0.0; % d5
pkin_gen(13) = 0.0; % d6
pkin_gen(14) = 0.0; % theta3
