% Inverse Dynamik für komplette Trajektorie für
% S6PRRRRR1
%
% Eingabe:
% RV_Traj [NTx109]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see S6PRRRRR1_invdynJ_fixb_regmin2vec.m
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR1_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx6]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S6PRRRRR1_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,109]), zeros(32,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 109), ...
  'S6PRRRRR1_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx109] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR1_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [32x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 6);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = S6PRRRRR1_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
