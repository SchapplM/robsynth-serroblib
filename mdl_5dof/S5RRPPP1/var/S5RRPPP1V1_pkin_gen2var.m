% Umwandlung der Kinematikparameter von S5RRPPP1 zu S5RRPPP1V1
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRPPP1
%   pkin_gen=[a2 a3 a4 a5 alpha3 d1 d2 theta3]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5RRPPP1V1
%   pkin_var=[a3 a4 a5 alpha3 d1 theta3]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRPPP1V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 6, 8];
pkin_var = pkin_gen(I_gv);
