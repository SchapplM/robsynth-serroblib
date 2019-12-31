% Umwandlung der Kinematikparameter von S6RRRPPP1 zu S6RRRPPP1V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRPPP1
%   pkin_gen=[a2 a3 a4 a5 a6 alpha4 d1 d2 d3 theta4]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRPPP1V1
%   pkin_var=[a4 a5 a6 alpha4 d1 theta4]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPPP1V1_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 5, 6, 7, 10];
pkin_var = pkin_gen(I_gv);
