% Umwandlung der Kinematikparameter von S6RRPRPR6 zu S6RRPRPR6V1
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPRPR6
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d4 d6 theta3]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRPRPR6V1
%   pkin_var=[a3 a4 a5 a6 d1 d4 d6 theta3]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRPR6V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 7, 9, 10, 11];
pkin_var = pkin_gen(I_gv);
