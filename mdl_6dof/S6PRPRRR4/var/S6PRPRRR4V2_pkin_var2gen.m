% Umwandlung der Kinematikparameter von S6PRPRRR4V2 zu S6PRPRRR4
% Eingabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6PRPRRR4V2
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRPRRR4
function pkin_gen = S6PRPRRR4V2_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 3, 6, 7, 8, 11, 12]) = pkin_var;
pkin_gen(4) = 0.0;
pkin_gen(5) = 0.0;
pkin_gen(9) = 0.0;
pkin_gen(10) = 0.0;
