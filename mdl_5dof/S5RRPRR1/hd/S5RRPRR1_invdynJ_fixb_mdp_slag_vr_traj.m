% Inverse Dynamik für komplette Trajektorie für
% S5RRPRR1
%
% Eingabe:
% RV_Traj [NTx74]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see S5RRPRR1_invdynJ_fixb_regmin2vec.m
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR1_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx5]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S5RRPRR1_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,74]), zeros(26,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 74), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx74] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [26x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 5);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = S5RRPRR1_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
