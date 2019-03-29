% Inverse Dynamik für komplette Trajektorie für
% S5RRRRR2
%
% Eingabe:
% RV_Traj [NTx89]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see S5RRRRR2_invdynJ_fixb_regmin2vec.m
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR2_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx5]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S5RRRRR2_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,89]), zeros(27,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 89), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx89] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [27x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 5);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = S5RRRRR2_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
