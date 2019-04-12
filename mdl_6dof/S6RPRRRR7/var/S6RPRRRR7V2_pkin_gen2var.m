% Umwandlung der Kinematikparameter von S6RPRRRR7 zu S6RPRRRR7V2
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RPRRRR7
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RPRRRR7V2
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRRR7V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8];
pkin_var = pkin_gen(I_gv);
