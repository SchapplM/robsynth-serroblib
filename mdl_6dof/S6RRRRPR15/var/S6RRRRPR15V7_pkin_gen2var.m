% Umwandlung der Kinematikparameter von S6RRRRPR15 zu S6RRRRPR15V7
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRPR15
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d6]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S6RRRRPR15V7
%   pkin_var=[a2 alpha2 d1 d2]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPR15V7_pkin_gen2var(pkin_gen)
I_gv = [1, 6, 8, 9];
pkin_var = pkin_gen(I_gv);
