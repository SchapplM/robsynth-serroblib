% Umwandlung der Kinematikparameter von S6RPRRRR11 zu S6RPRRRR11V2
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RPRRRR11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d3 d4 d5 d6 theta2]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RPRRRR11V2
%   pkin_var=[a2 a3 a6 alpha2 alpha3 d1 d3 d6 theta2]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRRR11V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 5, 6, 7, 8, 9, 12, 13];
pkin_var = pkin_gen(I_gv);
