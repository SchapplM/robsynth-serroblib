% Umwandlung der Kinematikparameter von S5RRPRR1 zu S5RRPRR1V1
% Eingabe:
% pkin_gen (4x1) double
%   Kinematikparameter (pkin) von S5RRPRR1
%   pkin_gen=[a3 a4 d4 d5]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRPRR1V1
%   pkin_var=[a3 a4 d4]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRR1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3];
pkin_var = pkin_gen(I_gv);
