% Umwandlung der Kinematikparameter von S6RRRRRR2V2 zu S6RRRRRR2
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRRR2V2
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRRRR2
function pkin_gen = S6RRRRRR2V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 3, 6, 7, 8, 9]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(10) = 0.0;
pkin_gen(11) = 0.0;
