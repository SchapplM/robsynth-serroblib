% Umwandlung der Kinematikparameter von S6RPPRRR5 zu S6RPPRRR5V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RPPRRR5
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RPPRRR5V2
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPPRRR5V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7];
pkin_var = pkin_gen(I_gv);
