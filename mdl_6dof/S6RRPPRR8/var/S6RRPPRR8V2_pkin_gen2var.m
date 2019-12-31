% Umwandlung der Kinematikparameter von S6RRPPRR8 zu S6RRPPRR8V2
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRPPRR8
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d5 d6 theta3]
% Ausgabe:
% pkin_var (8x1) double
%   Kinematikparameter (pkin) von S6RRPPRR8V2
%   pkin_var=[a3 a4 a5 a6 d1 d5 d6 theta3]
% I_gv (8x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPPRR8V2_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 6, 8, 9, 10];
pkin_var = pkin_gen(I_gv);
