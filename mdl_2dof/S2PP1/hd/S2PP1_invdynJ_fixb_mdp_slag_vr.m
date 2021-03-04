% Calculate vector of inverse dynamics joint torques for
% S2PP1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [3x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S2PP1_invdynJ_fixb_regmin2vec.m
% MDP [2x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2PP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-03-03 18:41
% Revision: 33b345ae0dd6ec4aa15499ab3d43edbbded0bea5 (2021-02-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S2PP1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1), zeros(2,1)}
assert(isreal(MDP) && all(size(MDP) == [2 1]), ...
  'S2PP1_invdynJ_fixb_mdp_slag_vr: MDP has to be [2x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2021-03-03 18:41:24
% EndTime: 2021-03-03 18:41:24
% DurationCPUTime: 0.02s
% Computational Cost: add. (1->1), mult. (3->3), div. (0->0), fcn. (3->3), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2); RV(3) * MDP(2);];
tauJ = t1;
