% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [93x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S6PRRPRR1_invdynJ_fixb_regmin2vec.m
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRR1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(93,1), zeros(29,1)}
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vr: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:29:06
% EndTime: 2021-01-16 03:29:06
% DurationCPUTime: 0.04s
% Computational Cost: add. (87->87), mult. (93->93), div. (0->0), fcn. (93->93), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(3) * MDP(3) + RV(5) * MDP(4) + RV(16) * MDP(10) + RV(19) * MDP(11) + RV(22) * MDP(12) + RV(26) * MDP(13) + RV(30) * MDP(14) + RV(34) * MDP(15) + RV(52) * MDP(21) + RV(57) * MDP(22) + RV(82) * MDP(28) + RV(88) * MDP(29); RV(27) * MDP(13) + RV(31) * MDP(14) + RV(35) * MDP(15) + RV(38) * MDP(16) + RV(41) * MDP(17) + RV(44) * MDP(18) + RV(47) * MDP(19) + RV(2) * MDP(2) + RV(4) * MDP(3) + RV(6) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(53) * MDP(21) + RV(58) * MDP(22) + RV(78) * MDP(27) + RV(83) * MDP(28) + RV(62) * MDP(23) + RV(66) * MDP(24) + RV(70) * MDP(25) + RV(74) * MDP(26) + RV(13) * MDP(8) + RV(17) * MDP(10) + RV(20) * MDP(11) + RV(23) * MDP(12) + RV(89) * MDP(29); RV(28) * MDP(13) + RV(32) * MDP(14) + RV(36) * MDP(15) + RV(39) * MDP(16) + RV(42) * MDP(17) + RV(45) * MDP(18) + RV(48) * MDP(19) + RV(8) * MDP(5) + RV(10) * MDP(6) + RV(50) * MDP(20) + RV(54) * MDP(21) + RV(59) * MDP(22) + RV(75) * MDP(26) + RV(79) * MDP(27) + RV(84) * MDP(28) + RV(63) * MDP(23) + RV(67) * MDP(24) + RV(71) * MDP(25) + RV(12) * MDP(7) + RV(14) * MDP(8) + RV(15) * MDP(9) + RV(18) * MDP(10) + RV(21) * MDP(11) + RV(24) * MDP(12) + RV(90) * MDP(29); RV(25) * MDP(12) + RV(29) * MDP(13) + RV(33) * MDP(14) + RV(37) * MDP(15) + RV(55) * MDP(21) + RV(60) * MDP(22) + RV(85) * MDP(28) + RV(91) * MDP(29); RV(40) * MDP(16) + RV(43) * MDP(17) + RV(46) * MDP(18) + RV(49) * MDP(19) + RV(51) * MDP(20) + RV(56) * MDP(21) + RV(61) * MDP(22) + RV(64) * MDP(23) + RV(68) * MDP(24) + RV(72) * MDP(25) + RV(76) * MDP(26) + RV(80) * MDP(27) + RV(86) * MDP(28) + RV(92) * MDP(29); RV(65) * MDP(23) + RV(69) * MDP(24) + RV(73) * MDP(25) + RV(77) * MDP(26) + RV(81) * MDP(27) + RV(87) * MDP(28) + RV(93) * MDP(29);];
tauJ = t1;
