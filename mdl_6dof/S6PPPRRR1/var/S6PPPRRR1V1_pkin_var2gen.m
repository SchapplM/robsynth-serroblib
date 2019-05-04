% Umwandlung der Kinematikparameter von S6PPPRRR1V1 zu S6PPPRRR1
% Eingabe:
% pkin_var (13x1) double
%   Kinematikparameter (pkin) von S6PPPRRR1V1
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6PPPRRR1
function pkin_gen = S6PPPRRR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14]) = pkin_var;
pkin_gen(11) = 0.0;
