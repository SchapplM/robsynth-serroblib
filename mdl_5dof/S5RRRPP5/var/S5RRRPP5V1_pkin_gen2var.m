% Umwandlung der Kinematikparameter von S5RRRPP5 zu S5RRRPP5V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S5RRRPP5
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3]
% Ausgabe:
% pkin_var (5x1) double
%   Kinematikparameter (pkin) von S5RRRPP5V1
%   pkin_var=[a3 a4 a5 d1 d3]
% I_gv (5x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRPP5V1_pkin_gen2var(pkin_gen)
I_gv = [2, 3, 4, 5, 7];
pkin_var = pkin_gen(I_gv);
