% Calculate vector of inverse dynamics joint torques for
% S5RPPRP6
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [58x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5RPPRP6_invdynJ_fixb_regmin2vec.m
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRP6_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(58,1), zeros(21,1)}
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRP6_invdynJ_fixb_mdp_slag_vr: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(8) + RV(16) * MDP(9) + RV(19) * MDP(10) + RV(22) * MDP(11) + RV(24) * MDP(12) + RV(26) * MDP(13) + RV(28) * MDP(14) + RV(31) * MDP(16) + RV(35) * MDP(17) + RV(39) * MDP(18) + RV(44) * MDP(19) + RV(49) * MDP(20) + RV(54) * MDP(21); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(14) * MDP(8) + RV(17) * MDP(9) + RV(20) * MDP(10) + RV(32) * MDP(16) + RV(36) * MDP(17) + RV(40) * MDP(18) + RV(45) * MDP(19) + RV(50) * MDP(20) + RV(55) * MDP(21); RV(12) * MDP(7) + RV(15) * MDP(8) + RV(18) * MDP(9) + RV(21) * MDP(10) + RV(33) * MDP(16) + RV(37) * MDP(17) + RV(41) * MDP(18) + RV(46) * MDP(19) + RV(51) * MDP(20) + RV(56) * MDP(21); RV(23) * MDP(11) + RV(25) * MDP(12) + RV(27) * MDP(13) + RV(29) * MDP(14) + RV(30) * MDP(15) + RV(34) * MDP(16) + RV(38) * MDP(17) + RV(42) * MDP(18) + RV(47) * MDP(19) + RV(52) * MDP(20) + RV(57) * MDP(21); RV(43) * MDP(18) + RV(48) * MDP(19) + RV(53) * MDP(20) + RV(58) * MDP(21);];
tauJ = t1;
