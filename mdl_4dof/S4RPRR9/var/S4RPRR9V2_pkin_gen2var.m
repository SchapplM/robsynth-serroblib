% Umwandlung der Kinematikparameter von S4RPRR9 zu S4RPRR9V2
% Eingabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RPRR9
%   pkin_gen=[a2 a3 a4 d1 d3 d4]
% Ausgabe:
% pkin_var (1x1) double
%   Kinematikparameter (pkin) von S4RPRR9V2
%   pkin_var=[d1]
% I_gv (1x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RPRR9V2_pkin_gen2var(pkin_gen)
I_gv = [4];
pkin_var = pkin_gen(I_gv);
