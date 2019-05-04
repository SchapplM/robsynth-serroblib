% Umwandlung der Kinematikparameter von S6PRPRRR7V1 zu S6PRPRRR7
% Eingabe:
% pkin_var (13x1) double
%   Kinematikparameter (pkin) von S6PRPRRR7V1
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6PRPRRR7
function pkin_gen = S6PRPRRR7V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14]) = pkin_var;
pkin_gen(12) = 0.0;
