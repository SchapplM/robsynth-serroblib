% Umwandlung der Kinematikparameter von S5RRRRR12 zu S5RRRRR12V1
% Eingabe:
% pkin_gen (11x1) double
%   Kinematikparameter (pkin) von S5RRRRR12
%   pkin_gen=[a2 a3 a4 a5 alpha2 alpha3 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S5RRRRR12V1
%   pkin_var=[a2 a3 alpha2 alpha3 d1 d2 d3]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRRR12V1_pkin_gen2var(pkin_gen)
I_gv = [1, 2, 5, 6, 7, 8, 9];
pkin_var = pkin_gen(I_gv);
