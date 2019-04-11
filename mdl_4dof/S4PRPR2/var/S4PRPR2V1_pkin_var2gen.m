% Umwandlung der Kinematikparameter von S4PRPR2V1 zu S4PRPR2
% Eingabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S4PRPR2V1
% Ausgabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4PRPR2
function pkin_gen = S4PRPR2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(6,1);
pkin_gen([1, 2, 3, 4, 6]) = pkin_var;
pkin_gen(5) = 0.0;
