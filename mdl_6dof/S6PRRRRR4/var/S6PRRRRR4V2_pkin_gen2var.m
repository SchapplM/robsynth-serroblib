% Umwandlung der Kinematikparameter von S6PRRRRR4 zu S6PRRRRR4V2
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6PRRRRR4
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRRRRR4V2
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRRRRR4V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8, 9, 10, 13];
pkin_var = pkin_gen(I_gv);
