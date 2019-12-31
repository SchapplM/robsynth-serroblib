% Umwandlung der Kinematikparameter von S6PPRRRR1 zu S6PPRRRR1V2
% Eingabe:
% pkin_gen (13x1) double
%   Kinematikparameter (pkin) von S6PPRRRR1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d3 d4 d5 d6 theta1 theta2]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PPRRRR1V2
%   pkin_var=[a2 a3 a5 alpha2 alpha3 d3 d5 theta1 theta2]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PPRRRR1V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 6, 7, 8, 10, 12, 13];
pkin_var = pkin_gen(I_gv);
