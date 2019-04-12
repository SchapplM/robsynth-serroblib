% Umwandlung der Kinematikparameter von S6PPRRRR1V2 zu S6PPRRRR1
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PPRRRR1V2
% Ausgabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6PPRRRR1
function pkin_gen = S6PPRRRR1V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(13,1);
pkin_gen([1, 2, 3, 6, 7, 8, 9, 12, 13]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(10) = 0.0;
pkin_gen(11) = 0.0;
