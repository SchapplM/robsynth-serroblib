% Umwandlung der Kinematikparameter von S6RRRRRR10V1 zu S6RRRRRR10
% Eingabe:
% pkin_var (13x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10V1
% Ausgabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10
function pkin_gen = S6RRRRRR10V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(14,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]) = pkin_var;
pkin_gen(14) = 0.0;
