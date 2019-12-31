% Umwandlung der Kinematikparameter von S6RRPRRP7 zu S6RRPRRP7V1
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S6RRPRRP7
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d4 d5]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S6RRPRRP7V1
%   pkin_var=[a3 a4 a6 d1 d4]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRP7V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 5, 6, 8];
pkin_var = pkin_gen(I_gv);
