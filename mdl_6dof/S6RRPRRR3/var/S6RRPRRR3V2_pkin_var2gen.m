% Umwandlung der Kinematikparameter von S6RRPRRR3V2 zu S6RRPRRR3
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRPRRR3V2
% Ausgabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRRR3
function pkin_gen = S6RRPRRR3V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(11,1);
pkin_gen([1, 2, 3, 6, 7, 8, 11]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(9) = 0.0;
pkin_gen(10) = 0.0;
