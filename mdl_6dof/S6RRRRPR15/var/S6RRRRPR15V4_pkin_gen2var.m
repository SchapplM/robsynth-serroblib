% Umwandlung der Kinematikparameter von S6RRRRPR15 zu S6RRRRPR15V4
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRPR15
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d6]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPR15V4
%   pkin_var=[a3 a5 a6 alpha3 d1 d3 d6]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPR15V4_pkin_gen2var(pkin_gen)
I_gv = [2, 4, 5, 7, 8, 10, 12];
pkin_var = pkin_gen(I_gv);
