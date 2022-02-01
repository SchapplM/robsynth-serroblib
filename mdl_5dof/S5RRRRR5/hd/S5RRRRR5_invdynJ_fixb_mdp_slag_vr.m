% Calculate vector of inverse dynamics joint torques for
% S5RRRRR5
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [75x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5RRRRR5_invdynJ_fixb_regmin2vec.m
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR5_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(75,1), zeros(23,1)}
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vr: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:18
% EndTime: 2022-01-20 12:02:18
% DurationCPUTime: 0.04s
% Computational Cost: add. (70->70), mult. (75->75), div. (0->0), fcn. (75->75), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(8) + RV(16) * MDP(9) + RV(19) * MDP(10) + RV(23) * MDP(11) + RV(27) * MDP(12) + RV(31) * MDP(13) + RV(36) * MDP(15) + RV(40) * MDP(16) + RV(44) * MDP(17) + RV(49) * MDP(18) + RV(54) * MDP(19) + RV(59) * MDP(20) + RV(66) * MDP(22) + RV(71) * MDP(23); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(14) * MDP(8) + RV(17) * MDP(9) + RV(20) * MDP(10) + RV(24) * MDP(11) + RV(28) * MDP(12) + RV(32) * MDP(13) + RV(37) * MDP(15) + RV(41) * MDP(16) + RV(45) * MDP(17) + RV(50) * MDP(18) + RV(55) * MDP(19) + RV(60) * MDP(20) + RV(67) * MDP(22) + RV(72) * MDP(23); RV(12) * MDP(7) + RV(15) * MDP(8) + RV(18) * MDP(9) + RV(21) * MDP(10) + RV(25) * MDP(11) + RV(29) * MDP(12) + RV(33) * MDP(13) + RV(38) * MDP(15) + RV(42) * MDP(16) + RV(46) * MDP(17) + RV(51) * MDP(18) + RV(56) * MDP(19) + RV(61) * MDP(20) + RV(68) * MDP(22) + RV(73) * MDP(23); RV(22) * MDP(10) + RV(26) * MDP(11) + RV(30) * MDP(12) + RV(34) * MDP(13) + RV(35) * MDP(14) + RV(39) * MDP(15) + RV(43) * MDP(16) + RV(47) * MDP(17) + RV(52) * MDP(18) + RV(57) * MDP(19) + RV(62) * MDP(20) + RV(64) * MDP(21) + RV(69) * MDP(22) + RV(74) * MDP(23); RV(48) * MDP(17) + RV(53) * MDP(18) + RV(58) * MDP(19) + RV(63) * MDP(20) + RV(65) * MDP(21) + RV(70) * MDP(22) + RV(75) * MDP(23);];
tauJ = t1;
