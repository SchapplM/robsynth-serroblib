% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR3
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [87x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S6RRPRPR3_invdynJ_fixb_regmin2vec.m
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR3_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(87,1), zeros(28,1)}
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRPR3_invdynJ_fixb_mdp_slag_vr: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(72) * MDP(26) + RV(76) * MDP(27) + RV(82) * MDP(28) + RV(42) * MDP(19) + RV(46) * MDP(20) + RV(51) * MDP(21) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(9) + RV(15) * MDP(10) + RV(17) * MDP(11) + RV(20) * MDP(12) + RV(23) * MDP(13) + RV(56) * MDP(22) + RV(60) * MDP(23) + RV(64) * MDP(24) + RV(68) * MDP(25) + RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(26) * MDP(14) + RV(29) * MDP(15) + RV(32) * MDP(16) + RV(35) * MDP(17) + RV(38) * MDP(18); RV(73) * MDP(26) + RV(77) * MDP(27) + RV(83) * MDP(28) + RV(43) * MDP(19) + RV(47) * MDP(20) + RV(52) * MDP(21) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(12) * MDP(8) + RV(14) * MDP(9) + RV(16) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(12) + RV(24) * MDP(13) + RV(57) * MDP(22) + RV(61) * MDP(23) + RV(65) * MDP(24) + RV(69) * MDP(25) + RV(5) * MDP(4) + RV(7) * MDP(5) + RV(27) * MDP(14) + RV(30) * MDP(15) + RV(33) * MDP(16) + RV(36) * MDP(17) + RV(39) * MDP(18); RV(19) * MDP(11) + RV(22) * MDP(12) + RV(40) * MDP(18) + RV(44) * MDP(19) + RV(48) * MDP(20) + RV(53) * MDP(21) + RV(78) * MDP(27) + RV(84) * MDP(28); RV(25) * MDP(13) + RV(28) * MDP(14) + RV(31) * MDP(15) + RV(34) * MDP(16) + RV(37) * MDP(17) + RV(41) * MDP(18) + RV(45) * MDP(19) + RV(49) * MDP(20) + RV(54) * MDP(21) + RV(58) * MDP(22) + RV(62) * MDP(23) + RV(66) * MDP(24) + RV(70) * MDP(25) + RV(74) * MDP(26) + RV(79) * MDP(27) + RV(85) * MDP(28); RV(50) * MDP(20) + RV(55) * MDP(21) + RV(80) * MDP(27) + RV(86) * MDP(28); RV(59) * MDP(22) + RV(63) * MDP(23) + RV(67) * MDP(24) + RV(71) * MDP(25) + RV(75) * MDP(26) + RV(81) * MDP(27) + RV(87) * MDP(28);];
tauJ  = t1;
