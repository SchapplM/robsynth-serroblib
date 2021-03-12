% Umwandlung der Kinematikparameter von S6PRRRRR6 zu S6PRRRRR6V5
% Eingabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6PRRRRR6
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d2 d3 d4 d5 d6 theta1]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6PRRRRR6V5
%   pkin_var=[a2 a5 alpha2 d2 d5 theta1]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRRRRR6V5_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 6, 9, 12, 14];
pkin_var = pkin_gen(I_gv);
