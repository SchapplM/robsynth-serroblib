% Umwandlung der Kinematikparameter von S6RRPPRR2 zu S6RRPPRR2V2
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPPRR2
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d5 d6 theta3 theta4]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRPPRR2V2
%   pkin_var=[a3 a4 a5 a6 d1 d5 d6 theta3 theta4]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPPRR2V2_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 6, 8, 9, 10, 11];
pkin_var = pkin_gen(I_gv);
