% Umwandlung der Kinematikparameter von S6RPRPPR7 zu S6RPRPPR7V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPRPPR7
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RPRPPR7V1
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRPPR7V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 9];
pkin_var = pkin_gen(I_gv);
