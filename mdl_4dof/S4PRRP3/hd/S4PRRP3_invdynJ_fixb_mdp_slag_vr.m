% Calculate vector of inverse dynamics joint torques for
% S4PRRP3
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [34x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S4PRRP3_invdynJ_fixb_regmin2vec.m
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRP3_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(34,1), zeros(15,1)}
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP3_invdynJ_fixb_mdp_slag_vr: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:09
% EndTime: 2021-01-14 22:27:09
% DurationCPUTime: 0.02s
% Computational Cost: add. (30->30), mult. (34->34), div. (0->0), fcn. (34->34), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(14) * MDP(10) + RV(17) * MDP(11) + RV(20) * MDP(12) + RV(24) * MDP(13) + RV(31) * MDP(15); RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(5) * MDP(5) + RV(7) * MDP(6) + RV(9) * MDP(7) + RV(11) * MDP(8) + RV(15) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(12) + RV(25) * MDP(13) + RV(28) * MDP(14) + RV(32) * MDP(15); RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(12) * MDP(8) + RV(13) * MDP(9) + RV(16) * MDP(10) + RV(19) * MDP(11) + RV(22) * MDP(12) + RV(26) * MDP(13) + RV(29) * MDP(14) + RV(33) * MDP(15); RV(23) * MDP(12) + RV(27) * MDP(13) + RV(30) * MDP(14) + RV(34) * MDP(15);];
tauJ = t1;
