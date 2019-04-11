% Umwandlung der Kinematikparameter von S4RPPR1V1 zu S4RPPR1
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4RPPR1V1
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RPPR1
function pkin_gen = S4RPPR1V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([1, 2, 3, 4, 6]) = pkin_var;
pkin_gen(5) = 0.0;
