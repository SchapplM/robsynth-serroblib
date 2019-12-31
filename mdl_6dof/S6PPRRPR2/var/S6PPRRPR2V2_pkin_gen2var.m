% Umwandlung der Kinematikparameter von S6PPRRPR2 zu S6PPRRPR2V2
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PPRRPR2
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d3 d4 d6 theta1 theta2]
% Ausgabe:
% pkin_var (10x1) double
%   Kinematikparameter (pkin) von S6PPRRPR2V2
%   pkin_var=[a2 a3 a5 a6 alpha2 alpha3 d3 d6 theta1 theta2]
% I_gv (10x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PPRRPR2V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 5, 6, 7, 8, 10, 11, 12];
pkin_var = pkin_gen(I_gv);
