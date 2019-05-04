% Umwandlung der Kinematikparameter von S6RPPRRR5V2 zu S6RPPRRR5
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RPPRRR5V2
% Ausgabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPPRRR5
function pkin_gen = S6RPPRRR5V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(9,1);
pkin_gen([1, 2, 3, 6, 7]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(8) = 0.0;
pkin_gen(9) = 0.0;
