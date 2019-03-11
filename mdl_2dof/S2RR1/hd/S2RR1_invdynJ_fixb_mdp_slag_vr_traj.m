% Inverse Dynamik für komplette Trajektorie für
% S2RR1
%
% Eingabe:
% RV_Traj [NTx16]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see S2RR1_invdynJ_fixb_regmin2vec.m
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR1_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx2]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S2RR1_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,16]), zeros(10,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 16), ...
  'S2RR1_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx16] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S2RR1_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [10x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 2);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = S2RR1_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
