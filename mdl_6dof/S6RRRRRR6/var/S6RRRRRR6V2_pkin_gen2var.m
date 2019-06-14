% Umwandlung der Kinematikparameter von S6RRRRRR6 zu S6RRRRRR6V2
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRRR6
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRRRR6V2
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRR6V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8, 9, 10];
pkin_var = pkin_gen(I_gv);