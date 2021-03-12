% Umwandlung der Kinematikparameter von S6PRRRRR5 zu S6PRRRRR5V3
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6PRRRRR5
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d2 d3 d4 d5 d6 theta1]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6PRRRRR5V3
%   pkin_var=[a2 a5 a6 alpha2 d2 d5 d6 theta1]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRRRRR5V3_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 8, 11, 12, 13];
pkin_var = pkin_gen(I_gv);
