% Umwandlung der Kinematikparameter von S6RPRRRR10 zu S6RPRRRR10V5
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RPRRRR10
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d3 d4 d5 d6 theta2]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RPRRRR10V5
%   pkin_var=[a5 alpha2 alpha3 d1 d5 theta2]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRRR10V5_pkin_gen2var(pkin_gen)
I_gv = [4, 6, 7, 8, 11, 13];
pkin_var = pkin_gen(I_gv);
