% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR5
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [100x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S6RPRRPR5_invdynJ_fixb_regmin2vec.m
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR5_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(100,1), zeros(32,1)}
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR5_invdynJ_fixb_mdp_slag_vr: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(89) * MDP(31) + RV(95) * MDP(32) + RV(85) * MDP(30) + RV(81) * MDP(29) + RV(77) * MDP(28) + RV(73) * MDP(27) + RV(69) * MDP(26) + RV(64) * MDP(25) + RV(49) * MDP(22) + RV(54) * MDP(23) + RV(59) * MDP(24) + RV(33) * MDP(17) + RV(36) * MDP(18) + RV(41) * MDP(20) + RV(45) * MDP(21) + RV(14) * MDP(9) + RV(16) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(13) + RV(24) * MDP(14) + RV(27) * MDP(15) + RV(30) * MDP(16) + RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(12) * MDP(8); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(22) * MDP(13) + RV(25) * MDP(14) + RV(42) * MDP(20) + RV(46) * MDP(21) + RV(50) * MDP(22) + RV(55) * MDP(23) + RV(60) * MDP(24) + RV(65) * MDP(25) + RV(90) * MDP(31) + RV(96) * MDP(32); RV(91) * MDP(31) + RV(97) * MDP(32) + RV(86) * MDP(30) + RV(78) * MDP(28) + RV(82) * MDP(29) + RV(74) * MDP(27) + RV(70) * MDP(26) + RV(66) * MDP(25) + RV(61) * MDP(24) + RV(47) * MDP(21) + RV(51) * MDP(22) + RV(56) * MDP(23) + RV(31) * MDP(16) + RV(34) * MDP(17) + RV(37) * MDP(18) + RV(39) * MDP(19) + RV(43) * MDP(20) + RV(15) * MDP(9) + RV(17) * MDP(10) + RV(19) * MDP(11) + RV(20) * MDP(12) + RV(23) * MDP(13) + RV(26) * MDP(14) + RV(28) * MDP(15) + RV(13) * MDP(8); RV(29) * MDP(15) + RV(32) * MDP(16) + RV(35) * MDP(17) + RV(38) * MDP(18) + RV(40) * MDP(19) + RV(44) * MDP(20) + RV(48) * MDP(21) + RV(52) * MDP(22) + RV(57) * MDP(23) + RV(62) * MDP(24) + RV(67) * MDP(25) + RV(71) * MDP(26) + RV(75) * MDP(27) + RV(79) * MDP(28) + RV(83) * MDP(29) + RV(87) * MDP(30) + RV(92) * MDP(31) + RV(98) * MDP(32); RV(53) * MDP(22) + RV(58) * MDP(23) + RV(63) * MDP(24) + RV(68) * MDP(25) + RV(93) * MDP(31) + RV(99) * MDP(32); RV(72) * MDP(26) + RV(76) * MDP(27) + RV(80) * MDP(28) + RV(84) * MDP(29) + RV(88) * MDP(30) + RV(94) * MDP(31) + RV(100) * MDP(32);];
tauJ  = t1;
