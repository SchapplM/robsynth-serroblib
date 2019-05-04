% Umwandlung der Kinematikparameter von S6RPRRPR1V1 zu S6RPRRPR1
% Eingabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6RPRRPR1V1
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RPRRPR1
function pkin_gen = S6RPRRPR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 10, 11]) = pkin_var;
pkin_gen(9) = 0.0;
