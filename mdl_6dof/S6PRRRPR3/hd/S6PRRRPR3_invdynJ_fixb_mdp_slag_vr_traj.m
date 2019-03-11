% Inverse Dynamik für komplette Trajektorie für
% S6PRRRPR3
%
% Eingabe:
% RV_Traj [NTx95]
%   time series of regressor matrices as vectors
%   Number of time steps (NT) in rows
%   see S6PRRRPR3_invdynJ_fixb_regmin2vec.m
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR3_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx6]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S6PRRRPR3_invdynJ_fixb_mdp_slag_vr_traj(RV_Traj, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,95]), zeros(29,1)}
assert(isreal(RV_Traj) && all(size(RV_Traj,2) == 95), ...
  'S6PRRRPR3_invdynJ_fixb_mdp_slag_vr_traj: RV_Traj needs to be [NTx95] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR3_invdynJ_fixb_mdp_slag_vr_traj: Dynamics parameter vector MDP has to be [29x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(RV_Traj,1), 6);
for ii = 1:size(RV_Traj,1)
  TAU(ii,:) = S6PRRRPR3_invdynJ_fixb_mdp_slag_vr(RV_Traj(ii,:), MDP);
end
