% Umwandlung der Kinematikparameter von S6PRRPRP1 zu S6PRRPRP1V1
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6PRRPRP1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d2 d3 d5 theta1 theta4]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6PRRPRP1V1
%   pkin_var=[a2 a4 a5 a6 alpha2 d2 d5 theta1 theta4]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6PRRPRP1V1_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 5, 6, 7, 9, 10, 11];
pkin_var = pkin_gen(I_gv);
