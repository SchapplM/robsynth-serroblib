% Umwandlung der Kinematikparameter von S5RRRPR5 zu S5RRRPR5V2
% Eingabe:
% pkin_gen (9x1) double
%   Kinematikparameter (pkin) von S5RRRPR5
%   pkin_gen=[a2 a3 a4 a5 d1 d2 d3 d5 theta4]
% Ausgabe:
% pkin_var (4x1) double
%   Kinematikparameter (pkin) von S5RRRPR5V2
%   pkin_var=[a3 d1 d3 theta4]
% I_gv (4x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S5RRRPR5V2_pkin_gen2var(pkin_gen)
I_gv = [2, 5, 7, 9];
pkin_var = pkin_gen(I_gv);
