% Umwandlung der Kinematikparameter von S6RRRPPR5 zu S6RRRPPR5V3
% Eingabe:
% pkin_gen (12x1) double
%   Kinematikparameter (pkin) von S6RRRPPR5
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d6 theta4 theta5]
% Ausgabe:
% pkin_var (7x1) double
%   Kinematikparameter (pkin) von S6RRRPPR5V3
%   pkin_var=[a4 a5 a6 d1 d6 theta4 theta5]
% I_gv (7x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRPPR5V3_pkin_gen2var(pkin_gen)
I_gv = [3, 4, 5, 7, 10, 11, 12];
pkin_var = pkin_gen(I_gv);
