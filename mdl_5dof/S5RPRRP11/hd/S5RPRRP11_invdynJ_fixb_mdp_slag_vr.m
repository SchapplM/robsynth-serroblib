% Calculate vector of inverse dynamics joint torques for
% S5RPRRP11
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [69x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5RPRRP11_invdynJ_fixb_regmin2vec.m
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP11_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(69,1), zeros(25,1)}
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRP11_invdynJ_fixb_mdp_slag_vr: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(12) * MDP(8) + RV(14) * MDP(9) + RV(16) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(13) + RV(24) * MDP(14) + RV(27) * MDP(15) + RV(30) * MDP(16) + RV(33) * MDP(17) + RV(36) * MDP(18) + RV(39) * MDP(19) + RV(42) * MDP(20) + RV(46) * MDP(21) + RV(50) * MDP(22) + RV(55) * MDP(23) + RV(60) * MDP(24) + RV(65) * MDP(25); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(22) * MDP(13) + RV(25) * MDP(14) + RV(43) * MDP(20) + RV(47) * MDP(21) + RV(51) * MDP(22) + RV(56) * MDP(23) + RV(61) * MDP(24) + RV(66) * MDP(25); RV(13) * MDP(8) + RV(15) * MDP(9) + RV(17) * MDP(10) + RV(19) * MDP(11) + RV(20) * MDP(12) + RV(23) * MDP(13) + RV(26) * MDP(14) + RV(28) * MDP(15) + RV(31) * MDP(16) + RV(34) * MDP(17) + RV(37) * MDP(18) + RV(40) * MDP(19) + RV(44) * MDP(20) + RV(48) * MDP(21) + RV(52) * MDP(22) + RV(57) * MDP(23) + RV(62) * MDP(24) + RV(67) * MDP(25); RV(29) * MDP(15) + RV(32) * MDP(16) + RV(35) * MDP(17) + RV(38) * MDP(18) + RV(41) * MDP(19) + RV(45) * MDP(20) + RV(49) * MDP(21) + RV(53) * MDP(22) + RV(58) * MDP(23) + RV(63) * MDP(24) + RV(68) * MDP(25); RV(54) * MDP(22) + RV(59) * MDP(23) + RV(64) * MDP(24) + RV(69) * MDP(25);];
tauJ = t1;
