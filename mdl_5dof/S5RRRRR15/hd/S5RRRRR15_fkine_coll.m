% Direkte kinematik für seriellen Roboter mit EE-Transformation für Kollisionsprüfung
%
% Eingabe:
% q
%   Gelenkwinkel des Roboters
% pkin
%   Kinematik-Parameter
% T_N_E
%   Transformationsmatrix EE-Segment-KS -> EE-KS
% I_EElink
%   Nummer des Segmentes, an dem der EE befestigt ist (0=Basis)
%
% Ausgabe:
% Tc_stack [(5+1+1)*3 x 4]
%   Gestapelte homogene Transformationsmatrizen für q (jew. ohne 0001-Zeile)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-05
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function Tc_stack = S5RRRRR15_fkine_coll(q, pkin, T_N_E, I_EElink)

%% Init
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(4,4),uint8(0)}

%% Berechnung, siehe S5RRRRR15_constr2.m
Tc_stack = NaN(3*(6+1),4);
[Tc_ges, Tc_stack(1:end-3,:)] = S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar(q, pkin);
T_0_E_q = Tc_ges(:,:,I_EElink+1) * T_N_E;
Tc_stack(end-2:end,:) = T_0_E_q(1:3,:);

