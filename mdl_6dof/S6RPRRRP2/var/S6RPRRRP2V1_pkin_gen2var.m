% Umwandlung der Kinematikparameter von S6RPRRRP2 zu S6RPRRRP2V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RPRRRP2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d3 d4 d5 theta2]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RPRRRP2V1
%   pkin_var=[a2 a3 a5 a6 d1 d3 d5 theta2]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RPRRRP2V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 5, 6, 7, 9, 10];
pkin_var = pkin_gen(I_gv);
