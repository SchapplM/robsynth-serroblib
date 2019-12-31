% Umwandlung der Kinematikparameter von S6RRRPRP5 zu S6RRRPRP5V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRPRP5
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d5 theta4]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRPRP5V1
%   pkin_var=[a4 a5 a6 d1 d5 theta4]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPRP5V1_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 5, 6, 9, 10];
pkin_var = pkin_gen(I_gv);
