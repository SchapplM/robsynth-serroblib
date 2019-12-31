% Umwandlung der Kinematikparameter von S6RRRRPP1 zu S6RRRRPP1V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRRPP1
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 theta5]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRPP1V1
%   pkin_var=[a3 a5 a6 d1 d3 theta5]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPP1V1_pkin_gen2var(pkin_gen)
I_gv = [2, 4, 5, 6, 8, 10];
pkin_var = pkin_gen(I_gv);
