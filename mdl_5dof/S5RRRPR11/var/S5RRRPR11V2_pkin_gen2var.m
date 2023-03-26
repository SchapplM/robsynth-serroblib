% Umwandlung der Kinematikparameter von S5RRRPR11 zu S5RRRPR11V2
% Eingabe:
% pkin_gen (8x1) double
%   Kinematikparameter (pkin) von S5RRRPR11
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d5]
% Ausgabe:
% pkin_var (1x1) double
%   Kinematikparameter (pkin) von S5RRRPR11V2
%   pkin_var=[d1]
% I_gv (1x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRPR11V2_pkin_gen2var(pkin_gen)
I_gv = [5];
pkin_var = pkin_gen(I_gv);
