% Umwandlung der Kinematikparameter von S4RRRR5 zu S4RRRR5V1
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RRRR5
%   pkin_gen=[a2 a3 a4 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRRR5V1
%   pkin_var=[a4 d1 d4]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRRR5V1_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 7];
pkin_var = pkin_gen(I_gv);
