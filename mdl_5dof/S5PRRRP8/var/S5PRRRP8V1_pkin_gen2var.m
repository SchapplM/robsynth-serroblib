% Umwandlung der Kinematikparameter von S5PRRRP8 zu S5PRRRP8V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5PRRRP8
%   pkin_gen=[a2 a3 a4 a5 alpha2 d2 d3 d4 theta1]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5PRRRP8V1
%   pkin_var=[a2 a5 alpha2 d2 theta1]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PRRRP8V1_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 9];
pkin_var = pkin_gen(I_gv);
