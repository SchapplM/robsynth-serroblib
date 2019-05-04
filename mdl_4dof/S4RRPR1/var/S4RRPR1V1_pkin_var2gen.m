% Umwandlung der Kinematikparameter von S4RRPR1V1 zu S4RRPR1
% Eingabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S4RRPR1V1
% Ausgabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RRPR1
function pkin_gen = S4RRPR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(7,1);
pkin_gen([1, 2, 3, 4, 5, 7]) = pkin_var;
pkin_gen(6) = 0.0;
