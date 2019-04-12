% Umwandlung der Kinematikparameter von S6RRRRRR8V2 zu S6RRRRRR8
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8V2
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRRR8
function pkin_gen = S6RRRRRR8V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 2, 3, 6, 7, 8, 9, 10, 11]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(12) = 0.0;
pkin_gen(13) = 0.0;
