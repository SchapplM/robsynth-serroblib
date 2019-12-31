% Umwandlung der Kinematikparameter von S6PRRRRR3 zu S6PRRRRR3V2
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRRRRR3
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d3 d4 d5 d6 theta1]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6PRRRRR3V2
%   pkin_var=[a2 a5 a6 alpha2 d2 d5 d6 theta1]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRRRRR3V2_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 7, 10, 11, 12];
pkin_var = pkin_gen(I_gv);
