% Umwandlung der Kinematikparameter von S6RRPPPR2V1 zu S6RRPPPR2
% Eingabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRPPPR2V1
% Ausgabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPPPR2
function pkin_gen = S6RRPPPR2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(10,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 9, 10]) = pkin_var;
pkin_gen(8) = 0.0;
