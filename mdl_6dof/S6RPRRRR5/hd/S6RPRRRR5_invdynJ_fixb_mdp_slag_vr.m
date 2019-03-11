% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR5
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [115x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S6RPRRRR5_invdynJ_fixb_regmin2vec.m
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR5_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(115,1), zeros(35,1)}
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR5_invdynJ_fixb_mdp_slag_vr: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(104) * MDP(34) + RV(110) * MDP(35) + RV(89) * MDP(31) + RV(94) * MDP(32) + RV(99) * MDP(33) + RV(74) * MDP(28) + RV(79) * MDP(29) + RV(84) * MDP(30) + RV(57) * MDP(24) + RV(61) * MDP(25) + RV(65) * MDP(26) + RV(69) * MDP(27) + RV(45) * MDP(21) + RV(49) * MDP(22) + RV(53) * MDP(23) + RV(30) * MDP(16) + RV(33) * MDP(17) + RV(36) * MDP(18) + RV(41) * MDP(20) + RV(12) * MDP(8) + RV(14) * MDP(9) + RV(16) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(13) + RV(24) * MDP(14) + RV(27) * MDP(15) + RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(22) * MDP(13) + RV(25) * MDP(14) + RV(42) * MDP(20) + RV(46) * MDP(21) + RV(70) * MDP(27) + RV(75) * MDP(28) + RV(105) * MDP(34) + RV(111) * MDP(35); RV(106) * MDP(34) + RV(112) * MDP(35) + RV(90) * MDP(31) + RV(95) * MDP(32) + RV(100) * MDP(33) + RV(76) * MDP(28) + RV(80) * MDP(29) + RV(85) * MDP(30) + RV(58) * MDP(24) + RV(62) * MDP(25) + RV(66) * MDP(26) + RV(71) * MDP(27) + RV(47) * MDP(21) + RV(50) * MDP(22) + RV(54) * MDP(23) + RV(28) * MDP(15) + RV(31) * MDP(16) + RV(34) * MDP(17) + RV(37) * MDP(18) + RV(39) * MDP(19) + RV(43) * MDP(20) + RV(13) * MDP(8) + RV(15) * MDP(9) + RV(17) * MDP(10) + RV(19) * MDP(11) + RV(20) * MDP(12) + RV(23) * MDP(13) + RV(26) * MDP(14); RV(29) * MDP(15) + RV(32) * MDP(16) + RV(35) * MDP(17) + RV(38) * MDP(18) + RV(40) * MDP(19) + RV(44) * MDP(20) + RV(48) * MDP(21) + RV(51) * MDP(22) + RV(55) * MDP(23) + RV(59) * MDP(24) + RV(63) * MDP(25) + RV(67) * MDP(26) + RV(72) * MDP(27) + RV(77) * MDP(28) + RV(81) * MDP(29) + RV(86) * MDP(30) + RV(91) * MDP(31) + RV(96) * MDP(32) + RV(101) * MDP(33) + RV(107) * MDP(34) + RV(113) * MDP(35); RV(52) * MDP(22) + RV(56) * MDP(23) + RV(60) * MDP(24) + RV(64) * MDP(25) + RV(68) * MDP(26) + RV(73) * MDP(27) + RV(78) * MDP(28) + RV(82) * MDP(29) + RV(87) * MDP(30) + RV(92) * MDP(31) + RV(97) * MDP(32) + RV(102) * MDP(33) + RV(108) * MDP(34) + RV(114) * MDP(35); RV(83) * MDP(29) + RV(88) * MDP(30) + RV(93) * MDP(31) + RV(98) * MDP(32) + RV(103) * MDP(33) + RV(109) * MDP(34) + RV(115) * MDP(35);];
tauJ  = t1;
