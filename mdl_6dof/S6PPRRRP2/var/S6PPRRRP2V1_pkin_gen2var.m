% Umwandlung der Kinematikparameter von S6PPRRRP2 zu S6PPRRRP2V1
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6PPRRRP2
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d3 d4 d5 theta1 theta2]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6PPRRRP2V1
%   pkin_var=[a2 a3 a6 alpha2 alpha3 d3 theta1 theta2]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PPRRRP2V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 5, 6, 7, 8, 11, 12];
pkin_var = pkin_gen(I_gv);
