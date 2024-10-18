% Calculate vector of inverse dynamics joint torques for
% S5RRRRR13
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [81x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5RRRRR13_invdynJ_fixb_regmin2vec.m
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR13_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(81,1), zeros(23,1)}
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR13_invdynJ_fixb_mdp_slag_vr: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:24
% EndTime: 2024-09-27 17:32:24
% DurationCPUTime: 0.00s
% Computational Cost: add. (76->76), mult. (81->81), div. (0->0), fcn. (81->81), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(8) + RV(16) * MDP(9) + RV(19) * MDP(10) + RV(23) * MDP(11) + RV(27) * MDP(12) + RV(31) * MDP(13) + RV(35) * MDP(14) + RV(39) * MDP(15) + RV(43) * MDP(16) + RV(47) * MDP(17) + RV(52) * MDP(18) + RV(57) * MDP(19) + RV(62) * MDP(20) + RV(67) * MDP(21) + RV(72) * MDP(22) + RV(77) * MDP(23); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(14) * MDP(8) + RV(17) * MDP(9) + RV(20) * MDP(10) + RV(24) * MDP(11) + RV(28) * MDP(12) + RV(32) * MDP(13) + RV(36) * MDP(14) + RV(40) * MDP(15) + RV(44) * MDP(16) + RV(48) * MDP(17) + RV(53) * MDP(18) + RV(58) * MDP(19) + RV(63) * MDP(20) + RV(68) * MDP(21) + RV(73) * MDP(22) + RV(78) * MDP(23); RV(12) * MDP(7) + RV(15) * MDP(8) + RV(18) * MDP(9) + RV(21) * MDP(10) + RV(25) * MDP(11) + RV(29) * MDP(12) + RV(33) * MDP(13) + RV(37) * MDP(14) + RV(41) * MDP(15) + RV(45) * MDP(16) + RV(49) * MDP(17) + RV(54) * MDP(18) + RV(59) * MDP(19) + RV(64) * MDP(20) + RV(69) * MDP(21) + RV(74) * MDP(22) + RV(79) * MDP(23); RV(22) * MDP(10) + RV(26) * MDP(11) + RV(30) * MDP(12) + RV(34) * MDP(13) + RV(38) * MDP(14) + RV(42) * MDP(15) + RV(46) * MDP(16) + RV(50) * MDP(17) + RV(55) * MDP(18) + RV(60) * MDP(19) + RV(65) * MDP(20) + RV(70) * MDP(21) + RV(75) * MDP(22) + RV(80) * MDP(23); RV(51) * MDP(17) + RV(56) * MDP(18) + RV(61) * MDP(19) + RV(66) * MDP(20) + RV(71) * MDP(21) + RV(76) * MDP(22) + RV(81) * MDP(23);];
tauJ = t1;
