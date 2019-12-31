% Umwandlung der Kinematikparameter von S5RRRRP10 zu S5RRRRP10V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRRP10
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRRRP10V1
%   pkin_var=[a2 a5 alpha2 d1 d2]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRP10V1_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 7];
pkin_var = pkin_gen(I_gv);
