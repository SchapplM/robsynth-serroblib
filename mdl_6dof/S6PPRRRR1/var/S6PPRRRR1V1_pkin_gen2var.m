% Umwandlung der Kinematikparameter von S6PPRRRR1 zu S6PPRRRR1V1
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6PPRRRR1
% Ausgabe:
% pkin_var (12x1) double
%   Kinematikparameter (pkin) von S6PPRRRR1V1
% I_gv (12x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PPRRRR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13];
pkin_var = pkin_gen(I_gv);
