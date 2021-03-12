% Umwandlung der Kinematikparameter von S6PRRRPR8 zu S6PRRRPR8V3
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PRRRPR8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d2 d3 d4 d6 theta1]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6PRRRPR8V3
%   pkin_var=[a2 a5 a6 alpha2 d2 d6 theta1]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRRRPR8V3_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 8, 11, 12];
pkin_var = pkin_gen(I_gv);
