% Kinematische Zwangsbedingungen zwischen Ist- und Soll-Konfiguration
% Die Zwangsbedingungen geben die Abweichung zwischen einer Soll-Pose in
% EE-Koordinaten und der Ist-Pose aus gegebenen Gelenk-Koordinaten an.
% Vollständige Rotations- und Translationskomponenten
% Variante 2:
% * Vektor vom Basis- zum EE-KS (kein Unterschied zu SKM-Variante 1)
% * Absolute Rotation ausgedrückt in XYZ-Euler-Winkeln (entspricht PKM
%   Variante 2)
% * Rotationsfehler ausgedrückt in Euler-Winkeln (um raumfeste Achsen), je
%   nach Eingabeargument `reci` (entspricht teilweise PKM-Variante 2)
%
% Eingabe:
% q
%   Gelenkwinkel des Roboters
% xE
%   Endeffektorpose des Roboters bezüglich des Basis-KS
% pkin
%   Kinematik-Parameter
% T_N_E
%   Transformationsmatrix EE-Segment-KS -> EE-KS
% phiconv_W_E
%   Winkelkonvention der Euler-Winkel Welt->End-Effektor. Siehe eul2r.m
% I_EElink
%   Nummer des Segmentes, an dem der EE befestigt ist (0=Basis)
% reci
%   true: Nehme reziproke Euler-Winkel für Orientierungsfehler (z.B.
%   ZYX-Orientierungsfehler für XYZ-Absolutorientierung)
%   false: Gleiche Euler-Winkel für Fehler und Absolut [Standard]
%
% Ausgabe:
% Phi
%   Kinematische Zwangsbedingungen des Roboters: Maß für den Positions- und
%   Orientierungsfehler zwischen Ist-Pose aus gegebenen Gelenkwinkeln q und
%   Soll-Pose aus gegebenen EE-Koordinaten x

% Quellen:
% [A] Aufzeichnungen Schappler vom 15.06.2018
% [B] Aufzeichnungen Schappler vom 22.06.2018
% [C] Aufzeichnungen Schappler vom 27.07.2018

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-28 13:06
% Revision: 2bf3b907e1213de0593c9d1d0a7eb98ef6ddbfca (2019-02-28)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi = S6PRRRRP4_constr2(q, xE, pkin, T_N_E, phiconv_W_E, I_EElink, reci)

%% Init
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(4,4),uint8(2),uint8(0),true}

%% Allgemein
Tc_ges = S6PRRRRP4_fkine_fixb_rotmat_mdh_sym_varpar(q, pkin);
T_0_E_q = Tc_ges(:,:,I_EElink+1) * T_N_E;
% T_0_E_q = Rob.fkineEE(q);

%% Translatorischer Teil

r_0_E_x = xE(1:3);

% Direkte Kinematik der Beinkette
r_0_E_q = T_0_E_q(1:3,4);

% Gl. (A.23, B.22)
Phix = r_0_E_q - r_0_E_x;

%% Rotatorischer Teil
R_0_E_x = eul2r(xE(4:6), phiconv_W_E);
if reci
  [~,phiconv_delta] = euler_angle_properties(phiconv_W_E);
else
  phiconv_delta = phiconv_W_E;
end

R_0_E_q = T_0_E_q(1:3,1:3);

% Gl. (C.47)
R_Ex_Eq = R_0_E_x' * R_0_E_q;

% Differenz-Rotation mit z.B. mit ZYX-Euler-Winkel
% Gl. (C.30) (dort Vertauschung der Reihenfolge, hier nicht)
phiR = r2eul(R_Ex_Eq, phiconv_delta);

Phi = [Phix; phiR];