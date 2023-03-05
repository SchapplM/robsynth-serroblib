% Umwandlung der Kinematikparameter von S4RRRP2 zu S4RRRP2V2
% Eingabe:
% pkin_gen (6x1) double
%   Kinematikparameter (pkin) von S4RRRP2
%   pkin_gen=[a2 a3 a4 d1 d2 d3]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S4RRRP2V2
%   pkin_var=[a2 d1 d2]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S4RRRP2V2_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5];
pkin_var = pkin_gen(I_gv);
