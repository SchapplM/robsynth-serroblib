% Umwandlung der Kinematikparameter von S6RPPRRR6 zu S6RPPRRR6V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPPRRR6
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RPPRRR6V2
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPPRRR6V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7];
pkin_var = pkin_gen(I_gv);
