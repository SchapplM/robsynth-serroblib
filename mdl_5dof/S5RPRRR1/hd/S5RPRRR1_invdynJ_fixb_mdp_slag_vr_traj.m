% Inverse Dynamik für komplette Trajektorie für
% S5RPRRR1
%
% Eingabe:
% RV_Traj [NTxNOTDEFINED]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see S5RPRRR1_invdynJ_fixb_regmin2vec.m
% MDP [NOTDEFINEDx1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR1_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx5]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-12 14:37
% Revision: aab8d7cd0cba739f5e0ec8d53b8419901d1154b0 (2019-06-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S5RPRRR1_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,NOTDEFINED]), zeros(NOTDEFINED,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == NOTDEFINED), ...
  'S5RPRRR1_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTxNOTDEFINED] (double)');
assert(isreal(MDP) && all(size(MDP) == [NOTDEFINED 1]), ...
  'S5RPRRR1_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [NOTDEFINEDx1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 5);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = S5RPRRR1_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
