% Umwandlung der Kinematikparameter von S6RRRRPR15 zu S6RRRRPR15V10
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRPR15
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d6]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRRRPR15V10
%   pkin_var=[a2 a3 a4 alpha2 alpha3 d1 d2 d3 d4]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPR15V10_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8, 9, 10, 11];
pkin_var = pkin_gen(I_gv);
