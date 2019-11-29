% Umwandlung der Kinematikparameter von S6RRPRRR14 zu S6RRPRRR14V3
% Eingabe:
% pkin_gen (14x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14
%   pkin_gen=[a2 a3 a4 a5 a6 alpha2 alpha3 alpha4 d1 d2 d4 d5 d6 theta3]
% Ausgabe:
% pkin_var (0x1) double
%   Kinematikparameter (pkin) von S6RRPRRR14V3
%   pkin_var=[]
% I_gv (0x1)
%   Vektor mit Indizes zur Selektion von Kinematikparametern
function [pkin_var, I_gv] = S6RRPRRR14V3_pkin_gen2var(pkin_gen)
I_gv = [];
pkin_var = pkin_gen(I_gv);
