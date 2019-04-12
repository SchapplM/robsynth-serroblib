% Umwandlung der Kinematikparameter von S6RPRRRR3 zu S6RPRRRR3V2
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RPRRRR3
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPRRRR3V2
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRRR3V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8, 11];
pkin_var = pkin_gen(I_gv);
