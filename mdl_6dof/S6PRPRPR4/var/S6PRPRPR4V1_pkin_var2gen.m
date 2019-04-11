% Umwandlung der Kinematikparameter von S6PRPRPR4V1 zu S6PRPRPR4
% Eingabe:
% pkin_var (11x1) double
%   Kinematikparameter (pkin) von S6PRPRPR4V1
% Ausgabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRPRPR4
function pkin_gen = S6PRPRPR4V1_pkin_var2gen(pkin_var)
pkin_gen = zeros(12,1);
pkin_gen([1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12]) = pkin_var;
pkin_gen(9) = 0.0;
