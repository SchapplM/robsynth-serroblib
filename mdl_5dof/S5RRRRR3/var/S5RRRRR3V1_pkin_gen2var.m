% Umwandlung der Kinematikparameter von S5RRRRR3 zu S5RRRRR3V1
% Eingabe:
% pkin_gen (5x1) double
%   Kinematikparameter (pkin) von S5RRRRR3
%   pkin_gen=[a3 a4 a5 d1 d4]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRRRR3V1
%   pkin_var=[a3 a5 d1]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRR3V1_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4];
pkin_var = pkin_gen(I_gv);
