% Inverse Dynamik für komplette Trajektorie für
% S7RRRRRRR1
%
% Eingabe:
% RV_Traj [NTx191]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see S7RRRRRRR1_invdynJ_fixb_regmin2vec.m
% MDP [45x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S7RRRRRRR1_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx7]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S7RRRRRRR1_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,191]), zeros(45,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 191), ...
  'S7RRRRRRR1_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx191] (double)');
assert(isreal(MDP) && all(size(MDP) == [45 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [45x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 7);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = S7RRRRRRR1_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
