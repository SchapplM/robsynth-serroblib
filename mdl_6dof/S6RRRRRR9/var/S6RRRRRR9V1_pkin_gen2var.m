% Umwandlung der Kinematikparameter von S6RRRRRR9 zu S6RRRRRR9V1
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRRR9
% Ausgabe:
% pkin_var (12x1) double
%   Kinematikparameter (pkin) von S6RRRRRR9V1
% I_gv (12x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRR9V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
pkin_var = pkin_gen(I_gv);
