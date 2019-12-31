% Umwandlung der Kinematikparameter von S6RRRPPR1 zu S6RRRPPR1V2
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRPPR1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d6 theta4 theta5]
% Ausgabe:
% pkin_var (9x1) double
%   Kinematikparameter (pkin) von S6RRRPPR1V2
%   pkin_var=[a3 a4 a5 a6 d1 d3 d6 theta4 theta5]
% I_gv (9x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPPR1V2_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 6, 8, 9, 10, 11];
pkin_var = pkin_gen(I_gv);
