% Umwandlung der Kinematikparameter von S6RPPRPR6V1 zu S6RPPRPR6
% Eingabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPPRPR6V1
% Ausgabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S6RPPRPR6
function pkin_gen = S6RPPRPR6V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(8,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7]) = pkin_var;
pkin_gen(8) = 0.0;
