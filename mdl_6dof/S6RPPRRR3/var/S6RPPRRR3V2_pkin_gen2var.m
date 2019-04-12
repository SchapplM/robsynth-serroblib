% Umwandlung der Kinematikparameter von S6RPPRRR3 zu S6RPPRRR3V2
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RPPRRR3
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RPPRRR3V2
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPPRRR3V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 10];
pkin_var = pkin_gen(I_gv);
