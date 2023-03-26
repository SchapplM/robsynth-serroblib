% Umwandlung der Kinematikparameter von S6RPRRRR11 zu S6RPRRRR11V4
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RPRRRR11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d3 d4 d5 d6 theta2]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RPRRRR11V4
%   pkin_var=[a6 alpha2 alpha3 d1 d6 theta2]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRRR11V4_pkin_gen2var(pkin_gen)
I_gv = [5, 6, 7, 8, 12, 13];
pkin_var = pkin_gen(I_gv);
