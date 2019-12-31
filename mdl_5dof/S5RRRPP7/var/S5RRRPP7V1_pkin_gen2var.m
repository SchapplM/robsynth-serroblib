% Umwandlung der Kinematikparameter von S5RRRPP7 zu S5RRRPP7V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S5RRRPP7
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S5RRRPP7V1
%   pkin_var=[a4 a5 d1]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRPP7V1_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 5];
pkin_var = pkin_gen(I_gv);
