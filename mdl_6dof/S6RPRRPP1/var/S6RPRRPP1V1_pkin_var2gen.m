% Umwandlung der Kinematikparameter von S6RPRRPP1V1 zu S6RPRRPP1
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RPRRPP1V1
%   pkin_var=[a2 a3 a5 a6 d1 d3 theta2 theta5]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RPRRPP1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 theta2 theta5]
%
% Siehe auch: S6RPRRPP1_structural_kinematic_parameters.m
function pkin_gen = S6RPRRPP1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 4, 5, 6, 7, 9, 10]) = pkin_var;

pkin_gen(3) = 0.0; % a4
pkin_gen(8) = 0.0; % d4
