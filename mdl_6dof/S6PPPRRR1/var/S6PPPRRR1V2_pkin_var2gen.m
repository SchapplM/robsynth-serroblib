% Umwandlung der Kinematikparameter von S6PPPRRR1V2 zu S6PPPRRR1
% Eingabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6PPPRRR1V2
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6PPPRRR1
function pkin_gen = S6PPPRRR1V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([1, 2, 3, 6, 7, 8, 9, 12, 13, 14]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(10) = 0.0;
pkin_gen(11) = 0.0;
