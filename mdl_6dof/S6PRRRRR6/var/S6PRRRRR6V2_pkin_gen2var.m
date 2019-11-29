% Umwandlung der Kinematikparameter von S6PRRRRR6 zu S6PRRRRR6V2
% Eingabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6PRRRRR6
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d2 d3 d4 d5 d6 theta1]
% Ausgabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6PRRRRR6V2
%   pkin_var=[a2 a3 a4 alpha2 alpha3 alpha4 d2 d3 d4 theta1]
% I_gv (10x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRRRRR6V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8, 9, 10, 11, 14];
pkin_var = pkin_gen(I_gv);
