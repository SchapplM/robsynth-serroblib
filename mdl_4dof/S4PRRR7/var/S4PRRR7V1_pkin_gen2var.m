% Umwandlung der Kinematikparameter von S4PRRR7 zu S4PRRR7V1
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S4PRRR7
%   pkin_gen=[a2 a3 a4 alpha2 d2 d3 d4 theta1]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S4PRRR7V1
%   pkin_var=[a2 alpha2 d2 theta1]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4PRRR7V1_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 8];
pkin_var = pkin_gen(I_gv);
