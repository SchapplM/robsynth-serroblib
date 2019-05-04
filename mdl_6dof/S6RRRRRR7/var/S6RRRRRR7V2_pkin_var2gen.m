% Umwandlung der Kinematikparameter von S6RRRRRR7V2 zu S6RRRRRR7
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRRRR7V2
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRRR7
function pkin_gen = S6RRRRRR7V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 3, 6, 7, 8, 9, 10]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(11) = 0.0;
pkin_gen(12) = 0.0;
