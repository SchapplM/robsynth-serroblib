% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR7
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [93x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S6PRPRRR7_invdynJ_fixb_regmin2vec.m
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRPRRR7_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(93,1), zeros(29,1)}
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR7_invdynJ_fixb_mdp_slag_vr: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1) + RV(3) * MDP(3) + RV(5) * MDP(4) + RV(7) * MDP(5) + RV(10) * MDP(6) + RV(13) * MDP(7) + RV(16) * MDP(8) + RV(29) * MDP(14) + RV(33) * MDP(15) + RV(52) * MDP(21) + RV(57) * MDP(22) + RV(82) * MDP(28) + RV(88) * MDP(29); RV(89) * MDP(29) + RV(37) * MDP(16) + RV(40) * MDP(17) + RV(43) * MDP(18) + RV(46) * MDP(19) + RV(49) * MDP(20) + RV(2) * MDP(2) + RV(4) * MDP(3) + RV(6) * MDP(4) + RV(70) * MDP(25) + RV(74) * MDP(26) + RV(78) * MDP(27) + RV(83) * MDP(28) + RV(53) * MDP(21) + RV(58) * MDP(22) + RV(62) * MDP(23) + RV(66) * MDP(24) + RV(8) * MDP(5) + RV(11) * MDP(6) + RV(14) * MDP(7) + RV(17) * MDP(8) + RV(19) * MDP(9) + RV(21) * MDP(10) + RV(23) * MDP(11) + RV(25) * MDP(12) + RV(27) * MDP(13) + RV(30) * MDP(14) + RV(34) * MDP(15); RV(9) * MDP(5) + RV(12) * MDP(6) + RV(15) * MDP(7) + RV(18) * MDP(8) + RV(31) * MDP(14) + RV(35) * MDP(15) + RV(54) * MDP(21) + RV(59) * MDP(22) + RV(84) * MDP(28) + RV(90) * MDP(29); RV(20) * MDP(9) + RV(22) * MDP(10) + RV(24) * MDP(11) + RV(26) * MDP(12) + RV(28) * MDP(13) + RV(32) * MDP(14) + RV(36) * MDP(15) + RV(38) * MDP(16) + RV(41) * MDP(17) + RV(44) * MDP(18) + RV(47) * MDP(19) + RV(50) * MDP(20) + RV(55) * MDP(21) + RV(60) * MDP(22) + RV(63) * MDP(23) + RV(67) * MDP(24) + RV(71) * MDP(25) + RV(75) * MDP(26) + RV(79) * MDP(27) + RV(85) * MDP(28) + RV(91) * MDP(29); RV(39) * MDP(16) + RV(42) * MDP(17) + RV(45) * MDP(18) + RV(48) * MDP(19) + RV(51) * MDP(20) + RV(56) * MDP(21) + RV(61) * MDP(22) + RV(64) * MDP(23) + RV(68) * MDP(24) + RV(72) * MDP(25) + RV(76) * MDP(26) + RV(80) * MDP(27) + RV(86) * MDP(28) + RV(92) * MDP(29); RV(65) * MDP(23) + RV(69) * MDP(24) + RV(73) * MDP(25) + RV(77) * MDP(26) + RV(81) * MDP(27) + RV(87) * MDP(28) + RV(93) * MDP(29);];
tauJ  = t1;
