% Umwandlung der Kinematikparameter von S6RRRRPR5 zu S6RRRRPR5V4
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRRPR5
%   pkin_gen=[a2 a3 a4 a5 a6 d1 d2 d3 d4 d6]
% Ausgabe:
% pkin_var (3x1) double
%   Kinematikparameter (pkin) von S6RRRRPR5V4
%   pkin_var=[a3 d1 d3]
% I_gv (3x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPR5V4_pkin_gen2var(pkin_gen)
I_gv = [2, 6, 8];
pkin_var = pkin_gen(I_gv);
