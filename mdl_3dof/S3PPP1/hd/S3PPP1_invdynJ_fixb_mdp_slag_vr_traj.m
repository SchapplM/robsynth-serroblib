% Inverse Dynamik für komplette Trajektorie für
% S3PPP1
%
% Eingabe:
% RV_Traj [NTx6]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see S3PPP1_invdynJ_fixb_regmin2vec.m
% MDP [3x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PPP1_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx3]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-17 09:48
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S3PPP1_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,6]), zeros(3,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 6), ...
  'S3PPP1_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx6] (double)');
assert(isreal(MDP) && all(size(MDP) == [3 1]), ...
  'S3PPP1_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [3x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 3);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = S3PPP1_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
