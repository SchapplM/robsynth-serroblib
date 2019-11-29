% Umwandlung der Kinematikparameter von S6RRRRPR12 zu S6RRRRPR12V1
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRPR12
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d6 theta5]
% Ausgabe:
% pkin_var (12x1) double
%   Kinematikparameter (pkin) von S6RRRRPR12V1
%   pkin_var=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 theta5]
% I_gv (12x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPR12V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13];
pkin_var = pkin_gen(I_gv);
