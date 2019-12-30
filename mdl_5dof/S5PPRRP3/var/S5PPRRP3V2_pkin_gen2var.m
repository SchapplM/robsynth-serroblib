% Umwandlung der Kinematikparameter von S5PPRRP3 zu S5PPRRP3V2
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5PPRRP3
%   pkin_gen=[a2 a3 a4 a5 d3 d4 theta1 theta2]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5PPRRP3V2
%   pkin_var=[a2 a3 a4 a5 d3 d4 theta1]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PPRRP3V2_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 4, 5, 6, 7];
pkin_var = pkin_gen(I_gv);
