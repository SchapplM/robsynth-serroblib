% Umwandlung der Kinematikparameter von S5RRRRP10 zu S5RRRRP10V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRRP10
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRRRP10V2
%   pkin_var=[a4 a5 d1 d4]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRP10V2_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 6, 9];
pkin_var = pkin_gen(I_gv);
