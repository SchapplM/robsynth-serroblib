% Umwandlung der Kinematikparameter von S6RRRRPR7 zu S6RRRRPR7V4
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRPR7
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d4 d6 theta5]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPR7V4
%   pkin_var=[a2 a4 alpha2 d1 d2 d4 theta5]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPR7V4_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 6, 7, 8, 10, 12];
pkin_var = pkin_gen(I_gv);
