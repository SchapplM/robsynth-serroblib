% Calculate vector of inverse dynamics joint torques for
% S5PRRPR3
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [59x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5PRRPR3_invdynJ_fixb_regmin2vec.m
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPR3_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(59,1), zeros(22,1)}
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vr: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:22
% EndTime: 2021-01-15 15:42:22
% DurationCPUTime: 0.06s
% Computational Cost: add. (54->54), mult. (59->59), div. (0->0), fcn. (59->59), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(14) * MDP(10) + RV(17) * MDP(11) + RV(20) * MDP(12) + RV(24) * MDP(13) + RV(28) * MDP(14) + RV(32) * MDP(15) + RV(50) * MDP(21) + RV(55) * MDP(22); RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(5) * MDP(5) + RV(7) * MDP(6) + RV(9) * MDP(7) + RV(11) * MDP(8) + RV(15) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(12) + RV(25) * MDP(13) + RV(29) * MDP(14) + RV(33) * MDP(15) + RV(36) * MDP(16) + RV(39) * MDP(17) + RV(42) * MDP(18) + RV(45) * MDP(19) + RV(51) * MDP(21) + RV(56) * MDP(22); RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(12) * MDP(8) + RV(13) * MDP(9) + RV(16) * MDP(10) + RV(19) * MDP(11) + RV(22) * MDP(12) + RV(26) * MDP(13) + RV(30) * MDP(14) + RV(34) * MDP(15) + RV(37) * MDP(16) + RV(40) * MDP(17) + RV(43) * MDP(18) + RV(46) * MDP(19) + RV(48) * MDP(20) + RV(52) * MDP(21) + RV(57) * MDP(22); RV(23) * MDP(12) + RV(27) * MDP(13) + RV(31) * MDP(14) + RV(35) * MDP(15) + RV(53) * MDP(21) + RV(58) * MDP(22); RV(38) * MDP(16) + RV(41) * MDP(17) + RV(44) * MDP(18) + RV(47) * MDP(19) + RV(49) * MDP(20) + RV(54) * MDP(21) + RV(59) * MDP(22);];
tauJ = t1;
