% Umwandlung der Kinematikparameter von S5RRPRR16 zu S5RRPRR16V3
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR16
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d4 d5]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRPRR16V3
%   pkin_var=[a2 alpha2 d1 d2]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRR16V3_pkin_gen2var(pkin_gen)
I_gv = [1, 5, 6, 7];
pkin_var = pkin_gen(I_gv);
