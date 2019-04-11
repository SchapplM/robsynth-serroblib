% Umwandlung der Kinematikparameter von S6PRPRRR4V1 zu S6PRPRRR4
% Eingabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6PRPRRR4V1
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRPRRR4
function pkin_gen = S6PRPRRR4V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12]) = pkin_var;
pkin_gen(10) = 0.0;
