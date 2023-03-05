% Umwandlung der Kinematikparameter von S6RRRRPR1 zu S6RRRRPR1V4
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S6RRRRPR1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d6 theta5]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRRRPR1V4
%   pkin_var=[a2 a3 a4 d1 d2 d3 d4 theta5]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPR1V4_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 3, 6, 7, 8, 9, 11];
pkin_var = pkin_gen(I_gv);
