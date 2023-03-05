% Umwandlung der Kinematikparameter von S6RRRRPR5V6 zu S6RRRRPR5
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPR5V6
%   pkin_var=[a2 a3 a4 d1 d2 d3 d4]
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRRPR5
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d6]
%
% Siehe auch: S6RRRRPR5_structural_kinematic_parameters.m
function pkin_gen = S6RRRRPR5V6_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 6, 7, 8, 9]) = pkin_var;

pkin_gen(4) = 0.0; % a5
pkin_gen(5) = 0.0; % a6
pkin_gen(10) = 0.0; % d6
