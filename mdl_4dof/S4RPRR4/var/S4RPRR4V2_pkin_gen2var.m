% Umwandlung der Kinematikparameter von S4RPRR4 zu S4RPRR4V2
% Eingabe:
% pkin_gen (7x1) double
%   Kinematikparameter (pkin) von S4RPRR4
%   pkin_gen=[a2 a3 a4 d1 d3 d4 theta2]
% Ausgabe:
% pkin_var (2x1) double
%   Kinematikparameter (pkin) von S4RPRR4V2
%   pkin_var=[d1 theta2]
% I_gv (2x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RPRR4V2_pkin_gen2var(pkin_gen)
I_gv = [4, 7];
pkin_var = pkin_gen(I_gv);
