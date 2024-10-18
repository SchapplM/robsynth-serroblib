% Calculate vector of inverse dynamics joint torques for
% S5RRRRR14
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [93x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5RRRRR14_invdynJ_fixb_regmin2vec.m
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR14_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(93,1), zeros(27,1)}
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR14_invdynJ_fixb_mdp_slag_vr: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:49
% EndTime: 2024-09-27 18:43:49
% DurationCPUTime: 0.05s
% Computational Cost: add. (88->88), mult. (93->93), div. (0->0), fcn. (93->93), ass. (0->1)
t1 = [RV(64) * MDP(22) + RV(69) * MDP(23) + RV(74) * MDP(24) + RV(79) * MDP(25) + RV(84) * MDP(26) + RV(39) * MDP(16) + RV(43) * MDP(17) + RV(47) * MDP(18) + RV(51) * MDP(19) + RV(55) * MDP(20) + RV(59) * MDP(21) + RV(10) * MDP(7) + RV(13) * MDP(8) + RV(16) * MDP(9) + RV(19) * MDP(10) + RV(22) * MDP(11) + RV(25) * MDP(12) + RV(28) * MDP(13) + RV(31) * MDP(14) + RV(35) * MDP(15) + RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(89) * MDP(27); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(14) * MDP(8) + RV(17) * MDP(9) + RV(20) * MDP(10) + RV(23) * MDP(11) + RV(26) * MDP(12) + RV(29) * MDP(13) + RV(32) * MDP(14) + RV(36) * MDP(15) + RV(40) * MDP(16) + RV(44) * MDP(17) + RV(48) * MDP(18) + RV(52) * MDP(19) + RV(56) * MDP(20) + RV(60) * MDP(21) + RV(65) * MDP(22) + RV(70) * MDP(23) + RV(75) * MDP(24) + RV(80) * MDP(25) + RV(85) * MDP(26) + RV(90) * MDP(27); RV(12) * MDP(7) + RV(15) * MDP(8) + RV(18) * MDP(9) + RV(21) * MDP(10) + RV(24) * MDP(11) + RV(27) * MDP(12) + RV(30) * MDP(13) + RV(33) * MDP(14) + RV(37) * MDP(15) + RV(41) * MDP(16) + RV(45) * MDP(17) + RV(49) * MDP(18) + RV(53) * MDP(19) + RV(57) * MDP(20) + RV(61) * MDP(21) + RV(66) * MDP(22) + RV(71) * MDP(23) + RV(76) * MDP(24) + RV(81) * MDP(25) + RV(86) * MDP(26) + RV(91) * MDP(27); RV(34) * MDP(14) + RV(38) * MDP(15) + RV(42) * MDP(16) + RV(46) * MDP(17) + RV(50) * MDP(18) + RV(54) * MDP(19) + RV(58) * MDP(20) + RV(62) * MDP(21) + RV(67) * MDP(22) + RV(72) * MDP(23) + RV(77) * MDP(24) + RV(82) * MDP(25) + RV(87) * MDP(26) + RV(92) * MDP(27); RV(63) * MDP(21) + RV(68) * MDP(22) + RV(73) * MDP(23) + RV(78) * MDP(24) + RV(83) * MDP(25) + RV(88) * MDP(26) + RV(93) * MDP(27);];
tauJ = t1;
