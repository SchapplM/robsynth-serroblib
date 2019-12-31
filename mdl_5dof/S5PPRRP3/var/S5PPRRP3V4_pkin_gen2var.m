% Umwandlung der Kinematikparameter von S5PPRRP3 zu S5PPRRP3V4
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5PPRRP3
%   pkin_gen=[a2 a3 a4 a5 d3 d4 theta1 theta2]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5PPRRP3V4
%   pkin_var=[a2 a3 a5 d3 theta1]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PPRRP3V4_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 4, 5, 7];
pkin_var = pkin_gen(I_gv);
