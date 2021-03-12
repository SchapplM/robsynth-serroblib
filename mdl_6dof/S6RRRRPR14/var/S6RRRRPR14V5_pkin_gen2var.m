% Umwandlung der Kinematikparameter von S6RRRRPR14 zu S6RRRRPR14V5
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6RRRRPR14
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d6 theta5]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRPR14V5
%   pkin_var=[a4 a5 a6 d1 d4 d6 theta5]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPR14V5_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 5, 8, 11, 12, 13];
pkin_var = pkin_gen(I_gv);
