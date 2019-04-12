% Umwandlung der Kinematikparameter von S6PRPRRR7 zu S6PRPRRR7V2
% Eingabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6PRPRRR7
% Ausgabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6PRPRRR7V2
% I_gv (10x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRPRRR7V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8, 9, 10, 13, 14];
pkin_var = pkin_gen(I_gv);
