% Umwandlung der Kinematikparameter von S5RRRRR11 zu S5RRRRR11V2
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S5RRRRR11
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRRRR11V2
%   pkin_var=[a4 a5 d1 d4 d5]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRR11V2_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 6, 9, 10];
pkin_var = pkin_gen(I_gv);
