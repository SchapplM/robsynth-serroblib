% Umwandlung der Kinematikparameter von S6RRPPRR11 zu S6RRPPRR11V3
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRPPRR11
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d5 d6 theta4]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRPPRR11V3
%   pkin_var=[a3 a4 a5 d1 d5 theta4]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPPRR11V3_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 7, 9, 11];
pkin_var = pkin_gen(I_gv);
