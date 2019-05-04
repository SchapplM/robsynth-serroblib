% Umwandlung der Kinematikparameter von S6RPPRRR4V2 zu S6RPPRRR4
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RPPRRR4V2
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RPPRRR4
function pkin_gen = S6RPPRRR4V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 6, 7, 10]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(8) = 0.0;
pkin_gen(9) = 0.0;
