% Calculate vector of inverse dynamics joint torques for
% S3PPR1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [10x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S3PPR1_invdynJ_fixb_regmin2vec.m
% MDP [5x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S3PPR1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1), zeros(5,1)}
assert(isreal(MDP) && all(size(MDP) == [5 1]), ...
  'S3PPR1_invdynJ_fixb_mdp_slag_vr: MDP has to be [5x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(5) * MDP(4) + RV(8) * MDP(5); RV(3) * MDP(2) + RV(6) * MDP(4) + RV(9) * MDP(5); RV(4) * MDP(3) + RV(7) * MDP(4) + RV(10) * MDP(5);];
tauJ  = t1;
