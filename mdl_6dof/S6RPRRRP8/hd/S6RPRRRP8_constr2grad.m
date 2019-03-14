% Ableitung der kinematischen Zwangsbedingungen nach den Gelenkwinkeln
% Die Zwangsbedingungen geben die Abweichung zwischen einer Soll-Pose in
% EE-Koordinaten und der Ist-Pose aus gegebenen Gelenk-Koordinaten an.
% Variante 2:
% * Vektor vom Basis- zum EE-KS (kein Unterschied zu SKM-Variante 1)
% * Absolute Rotation ausgedrückt in XYZ-Euler-Winkeln (entspricht PKM
%   Variante 2)
% * Rotationsfehler ausgedrückt in Euler-Winkeln (um raumfeste Achsen), je
%   nach Eingabeargument `reci` (entspricht teilweise PKM-Variante 2)
%
% Eingabe:
% q
%   Gelenkkoordinaten des Roboters
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
% Phi_dq [6xN]
%   Matrix mit Ableitungen der 6 Zwangsbedingungskomponenten (in den Zeilen)
%   nach den N Gelenkwinkeln (in den Spalten)

% Quellen:
% [A] Aufzeichnungen Schappler vom 21.06.2018
% [B] Aufzeichnungen Schappler vom 13.07.2018
% [C] Aufzeichnungen Schappler vom 27.07.2018
% [D] Aufzeichnungen Schappler vom 21.08.2018

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-28 14:05
% Revision: 2bf3b907e1213de0593c9d1d0a7eb98ef6ddbfca (2019-02-28)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-07
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Phi_dq = S6RPRRRP8_constr2grad(q, xE, pkin, T_N_E, phiconv_W_E, I_EElink, reci)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(4,4),uint8(2),uint8(0), true}

%% Translatorisch
% Bein-Jacobi
J0_i_trans = S6RPRRRP8_jacobia_transl_sym_varpar(q, I_EElink, T_N_E(1:3,4), pkin);
J_Ai_Bi = J0_i_trans; % Nur xyz-Koordinate in ZB.
dPhit_dq = J_Ai_Bi;

%% Rotatorisch
R_0_E_x = eul2r(xE(4:6), phiconv_W_E);
if reci
  [~,phiconv_delta] = euler_angle_properties(phiconv_W_E);
else
  phiconv_delta = phiconv_W_E;
end

% Kinematik, Definitionen
Tc_ges = S6RPRRRP8_fkine_fixb_rotmat_mdh_sym_varpar(q, pkin);
R_0_E_q = Tc_ges(1:3,1:3,I_EElink+1) * T_N_E(1:3,1:3);
R_Ex_Eq = R_0_E_x' * R_0_E_q;

% Term III aus Gl. (C.49): Ableitung der Rotationsmatrix R_0_E nach q
% Gl. D.7
b11=T_N_E(1,1);b12=T_N_E(1,2);b13=T_N_E(1,3);
b21=T_N_E(2,1);b22=T_N_E(2,2);b23=T_N_E(2,3);
b31=T_N_E(3,1);b32=T_N_E(3,2);b33=T_N_E(3,3);
dPidRb1 = [b11 0 0 b21 0 0 b31 0 0; 0 b11 0 0 b21 0 0 b31 0; 0 0 b11 0 0 b21 0 0 b31; b12 0 0 b22 0 0 b32 0 0; 0 b12 0 0 b22 0 0 b32 0; 0 0 b12 0 0 b22 0 0 b32; b13 0 0 b23 0 0 b33 0 0; 0 b13 0 0 b23 0 0 b33 0; 0 0 b13 0 0 b23 0 0 b33;];
dRb_0N_dq = S6RPRRRP8_jacobiR_rot_sym_varpar(q, I_EElink, pkin);
% Gl. D.7 für Term III in Gl. (C.49) einsetzen
dRb_0E_dq = dPidRb1 * dRb_0N_dq;

% Term II aus Gl. (C.49): Innere Ableitung des Matrix-Produktes
% aus ZB_diff_q_rmatvecprod_diff_rmatvec2_matlab
% (Matrix R_0_E_x wird transponiert in Vorlage eingesetzt)
a11=R_0_E_x(1,1);a12=R_0_E_x(2,1);a13=R_0_E_x(3,1);
a21=R_0_E_x(1,2);a22=R_0_E_x(2,2);a23=R_0_E_x(3,2);
a31=R_0_E_x(1,3);a32=R_0_E_x(2,3);a33=R_0_E_x(3,3);
dPi_dRb2 = [a11 a12 a13 0 0 0 0 0 0; a21 a22 a23 0 0 0 0 0 0; a31 a32 a33 0 0 0 0 0 0; 0 0 0 a11 a12 a13 0 0 0; 0 0 0 a21 a22 a23 0 0 0; 0 0 0 a31 a32 a33 0 0 0; 0 0 0 0 0 0 a11 a12 a13; 0 0 0 0 0 0 a21 a22 a23; 0 0 0 0 0 0 a31 a32 a33;];

% Term I aus Gl. (C.49): Ableitung der Euler-Winkel nach der Rot.-matrix
% (ZYX-Euler-Winkel des Orientierungsfehlers)
% Unabhängig vom Roboter (nur abhängig von Orientierungsdarstellung)
ddeltaR_dRb = eul_diff_rotmat(R_Ex_Eq,phiconv_delta);

% Gl. (C.49) (dort Vertauschung der Reihenfolge)
dPhir_dq = ddeltaR_dRb * dPi_dRb2 * dRb_0E_dq;

%% Ausgabe
% Vollständige Gradientenmatrix (3x Translation, 3x Rotation)
Phi_dq = [dPhit_dq; dPhir_dq];