% Inverse Dynamik für komplette Trajektorie für
% S6RRRRRR10V2
%
% Eingabe:
% RV_Traj [NTx141]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see S6RRRRRR10V2_invdynJ_fixb_regmin2vec.m
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR10V2_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx6]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S6RRRRRR10V2_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,141]), zeros(38,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 141), ...
  'S6RRRRRR10V2_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx141] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [38x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 6);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = S6RRRRRR10V2_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
