% Umwandlung der Kinematikparameter von S4PPPR2V1 zu S4PPPR2
% Eingabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4PPPR2V1
% Ausgabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S4PPPR2
function pkin_gen = S4PPPR2V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(5,1);
pkin_gen([1, 2, 3, 5]) = pkin_var;
pkin_gen(4) = 0.0;