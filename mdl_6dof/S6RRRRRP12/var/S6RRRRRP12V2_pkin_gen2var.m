% Umwandlung der Kinematikparameter von S6RRRRRP12 zu S6RRRRRP12V2
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRRRP12
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 d1 d2 d3 d4 d5]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRRRP12V2
%   pkin_var=[a2 a5 a6 alpha2 d1 d2 d5]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRRP12V2_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 8, 9, 12];
pkin_var = pkin_gen(I_gv);
