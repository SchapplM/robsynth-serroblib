% Umwandlung der Kinematikparameter von S6PRPRRR5V1 zu S6PRPRRR5
% Eingabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6PRPRRR5V1
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6PRPRRR5
function pkin_gen = S6PRPRRR5V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 9, 11]) = pkin_var;
pkin_gen(10) = 0.0;
