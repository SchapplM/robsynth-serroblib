% Umwandlung der Kinematikparameter von S4RRRR4 zu S4RRRR4V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RRRR4
%   pkin_gen=[a2 a3 a4 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRRR4V1
%   pkin_var=[a3 d1 d3]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRRR4V1_pkin_gen2var(pkin_gen)
I_gv = [2, 4, 6];
pkin_var = pkin_gen(I_gv);
