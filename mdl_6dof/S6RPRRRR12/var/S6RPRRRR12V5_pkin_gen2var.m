% Umwandlung der Kinematikparameter von S6RPRRRR12 zu S6RPRRRR12V5
% Eingabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RPRRRR12
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d1 d3 d4 d5 d6 theta2]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RPRRRR12V5
%   pkin_var=[a4 alpha2 alpha3 alpha4 d1 d4 theta2]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRRR12V5_pkin_gen2var(pkin_gen)
I_gv = [3, 6, 7, 8, 9, 11, 14];
pkin_var = pkin_gen(I_gv);
