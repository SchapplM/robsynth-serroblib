% Umwandlung der Kinematikparameter von S5RRPRR16 zu S5RRPRR16V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRPRR16
%   pkin_gen=[a2 a3 a4 a5 alpha2 d1 d2 d4 d5]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRPRR16V2
%   pkin_var=[a3 a4 d1 d4]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPRR16V2_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 6, 8];
pkin_var = pkin_gen(I_gv);
