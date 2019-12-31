% Umwandlung der Kinematikparameter von S6RRRRPP8 zu S6RRRRPP8V1
% Eingabe:
% pkin_gen (10x1) double
%   Kinematikparameter (pkin) von S6RRRRPP8
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 d1 d2 d3 d4]
% Ausgabe:
% pkin_var (6x1) double
%   Kinematikparameter (pkin) von S6RRRRPP8V1
%   pkin_var=[a2 a5 a6 alpha2 d1 d2]
% I_gv (6x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRRRPP8V1_pkin_gen2var(pkin_gen)
I_gv = [1, 4, 5, 6, 7, 8];
pkin_var = pkin_gen(I_gv);
