% Umwandlung der Kinematikparameter von S6PRPRRR3 zu S6PRPRRR3V2
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRPRRR3
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6PRPRRR3V2
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRPRRR3V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8, 11, 12];
pkin_var = pkin_gen(I_gv);
