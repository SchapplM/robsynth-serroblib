% Umwandlung der Kinematikparameter von S4RRPP4V1 zu S4RRPP4
% Eingabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRPP4V1
%   pkin_var=[a3 a4 d1]
% Ausgabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S4RRPP4
%   pkin_gen=[a2 a3 a4 d1 d2]
%
% Siehe auch: S4RRPP4_structural_kinematic_parameters.m
function pkin_gen = S4RRPP4V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(5,1);
pkin_gen([2, 3, 4]) = pkin_var;

pkin_gen(1) = 0.0; % a2
pkin_gen(5) = 0.0; % d2
