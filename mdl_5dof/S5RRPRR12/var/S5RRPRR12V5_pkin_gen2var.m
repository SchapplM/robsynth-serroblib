% Umwandlung der Kinematikparameter von S5RRPRR12 zu S5RRPRR12V5
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRPRR12
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d4 d5]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRPRR12V5
%   pkin_var=[a2 a5 d1 d2 d5]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRR12V5_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 8];
pkin_var = pkin_gen(I_gv);
