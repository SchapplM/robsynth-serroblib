% Umwandlung der Kinematikparameter von S6RRRRRR10 zu S6RRRRRR10V2
% Eingabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRRR10V2
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRR10V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 9, 12, 14];
pkin_var = pkin_gen(I_gv);
