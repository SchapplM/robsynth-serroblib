% Umwandlung der Kinematikparameter von S5PRRPP2 zu S5PRRPP2V1
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5PRRPP2
%   pkin_gen=[a2 a3 a4 a5 d2 d3 theta1 theta4]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S5PRRPP2V1
%   pkin_var=[a2 a4 a5 d2 theta1 theta4]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5PRRPP2V1_pkin_gen2var(pkin_gen)
I_gv = [1, 3, 4, 5, 7, 8];
pkin_var = pkin_gen(I_gv);
