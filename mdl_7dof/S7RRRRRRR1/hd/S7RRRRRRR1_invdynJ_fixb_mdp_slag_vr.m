% Calculate vector of inverse dynamics joint torques for
% S7RRRRRRR1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [191x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S7RRRRRRR1_invdynJ_fixb_regmin2vec.m
% MDP [45x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S7RRRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [7x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S7RRRRRRR1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(191,1), zeros(45,1)}
assert(isreal(MDP) && all(size(MDP) == [45 1]), ...
  'S7RRRRRRR1_invdynJ_fixb_mdp_slag_vr: MDP has to be [45x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(178) * MDP(44) + RV(185) * MDP(45) + RV(164) * MDP(42) + RV(171) * MDP(43) + RV(150) * MDP(40) + RV(157) * MDP(41) + RV(137) * MDP(38) + RV(143) * MDP(39) + RV(131) * MDP(37) + RV(125) * MDP(36) + RV(113) * MDP(34) + RV(119) * MDP(35) + RV(107) * MDP(33) + RV(96) * MDP(31) + RV(101) * MDP(32) + RV(86) * MDP(29) + RV(91) * MDP(30) + RV(81) * MDP(28) + RV(76) * MDP(27) + RV(71) * MDP(26) + RV(66) * MDP(25) + RV(58) * MDP(23) + RV(62) * MDP(24) + RV(46) * MDP(20) + RV(50) * MDP(21) + RV(54) * MDP(22) + RV(35) * MDP(17) + RV(38) * MDP(18) + RV(42) * MDP(19) + RV(23) * MDP(13) + RV(26) * MDP(14) + RV(29) * MDP(15) + RV(32) * MDP(16) + RV(10) * MDP(7) + RV(13) * MDP(9) + RV(15) * MDP(10) + RV(17) * MDP(11) + RV(20) * MDP(12) + RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6); RV(179) * MDP(44) + RV(186) * MDP(45) + RV(165) * MDP(42) + RV(172) * MDP(43) + RV(151) * MDP(40) + RV(158) * MDP(41) + RV(138) * MDP(38) + RV(144) * MDP(39) + RV(126) * MDP(36) + RV(132) * MDP(37) + RV(120) * MDP(35) + RV(114) * MDP(34) + RV(102) * MDP(32) + RV(108) * MDP(33) + RV(97) * MDP(31) + RV(87) * MDP(29) + RV(92) * MDP(30) + RV(82) * MDP(28) + RV(77) * MDP(27) + RV(72) * MDP(26) + RV(67) * MDP(25) + RV(59) * MDP(23) + RV(63) * MDP(24) + RV(47) * MDP(20) + RV(51) * MDP(21) + RV(55) * MDP(22) + RV(33) * MDP(16) + RV(36) * MDP(17) + RV(39) * MDP(18) + RV(43) * MDP(19) + RV(24) * MDP(13) + RV(27) * MDP(14) + RV(30) * MDP(15) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(12) * MDP(8) + RV(14) * MDP(9) + RV(16) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(12) + RV(5) * MDP(4) + RV(7) * MDP(5); RV(180) * MDP(44) + RV(187) * MDP(45) + RV(166) * MDP(42) + RV(173) * MDP(43) + RV(152) * MDP(40) + RV(159) * MDP(41) + RV(139) * MDP(38) + RV(145) * MDP(39) + RV(127) * MDP(36) + RV(133) * MDP(37) + RV(121) * MDP(35) + RV(115) * MDP(34) + RV(103) * MDP(32) + RV(109) * MDP(33) + RV(98) * MDP(31) + RV(88) * MDP(29) + RV(93) * MDP(30) + RV(83) * MDP(28) + RV(78) * MDP(27) + RV(73) * MDP(26) + RV(68) * MDP(25) + RV(56) * MDP(22) + RV(60) * MDP(23) + RV(64) * MDP(24) + RV(48) * MDP(20) + RV(52) * MDP(21) + RV(34) * MDP(16) + RV(37) * MDP(17) + RV(40) * MDP(18) + RV(44) * MDP(19) + RV(22) * MDP(12) + RV(25) * MDP(13) + RV(28) * MDP(14) + RV(31) * MDP(15) + RV(19) * MDP(11); RV(181) * MDP(44) + RV(188) * MDP(45) + RV(167) * MDP(42) + RV(174) * MDP(43) + RV(146) * MDP(39) + RV(153) * MDP(40) + RV(160) * MDP(41) + RV(140) * MDP(38) + RV(128) * MDP(36) + RV(134) * MDP(37) + RV(122) * MDP(35) + RV(116) * MDP(34) + RV(104) * MDP(32) + RV(110) * MDP(33) + RV(99) * MDP(31) + RV(89) * MDP(29) + RV(94) * MDP(30) + RV(84) * MDP(28) + RV(79) * MDP(27) + RV(74) * MDP(26) + RV(65) * MDP(24) + RV(69) * MDP(25) + RV(57) * MDP(22) + RV(61) * MDP(23) + RV(49) * MDP(20) + RV(53) * MDP(21) + RV(41) * MDP(18) + RV(45) * MDP(19); RV(70) * MDP(25) + RV(75) * MDP(26) + RV(80) * MDP(27) + RV(85) * MDP(28) + RV(90) * MDP(29) + RV(95) * MDP(30) + RV(100) * MDP(31) + RV(105) * MDP(32) + RV(111) * MDP(33) + RV(117) * MDP(34) + RV(123) * MDP(35) + RV(129) * MDP(36) + RV(135) * MDP(37) + RV(141) * MDP(38) + RV(147) * MDP(39) + RV(154) * MDP(40) + RV(161) * MDP(41) + RV(168) * MDP(42) + RV(175) * MDP(43) + RV(182) * MDP(44) + RV(189) * MDP(45); RV(106) * MDP(32) + RV(112) * MDP(33) + RV(118) * MDP(34) + RV(124) * MDP(35) + RV(130) * MDP(36) + RV(136) * MDP(37) + RV(142) * MDP(38) + RV(148) * MDP(39) + RV(155) * MDP(40) + RV(162) * MDP(41) + RV(169) * MDP(42) + RV(176) * MDP(43) + RV(183) * MDP(44) + RV(190) * MDP(45); RV(149) * MDP(39) + RV(156) * MDP(40) + RV(163) * MDP(41) + RV(170) * MDP(42) + RV(177) * MDP(43) + RV(184) * MDP(44) + RV(191) * MDP(45);];
tauJ  = t1;
