% Calculate vector of inverse dynamics joint torques for
% S3PPP1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [6x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S3PPP1_invdynJ_fixb_regmin2vec.m
% MDP [3x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-17 09:48
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S3PPP1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1), zeros(3,1)}
assert(isreal(MDP) && all(size(MDP) == [3 1]), ...
  'S3PPP1_invdynJ_fixb_mdp_slag_vr: MDP has to be [3x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(4) * MDP(3); RV(3) * MDP(2) + RV(5) * MDP(3); RV(6) * MDP(3);];
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(4) * MDP(3); RV(3) * MDP(2) + RV(5) * MDP(3); RV(6) * MDP(3);];
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(4) * MDP(3); RV(3) * MDP(2) + RV(5) * MDP(3); RV(6) * MDP(3);];
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(4) * MDP(3); RV(3) * MDP(2) + RV(5) * MDP(3); RV(6) * MDP(3);];
tauJ  = t1;
