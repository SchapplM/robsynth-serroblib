% Calculate vector of inverse dynamics joint torques for
% S5RRPPR1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [45x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5RRPPR1_invdynJ_fixb_regmin2vec.m
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(45,1), zeros(17,1)}
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vr: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:49
% EndTime: 2022-01-20 09:51:49
% DurationCPUTime: 0.03s
% Computational Cost: add. (40->40), mult. (45->45), div. (0->0), fcn. (45->45), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(8) + RV(16) * MDP(9) + RV(19) * MDP(10) + RV(23) * MDP(11) + RV(26) * MDP(12) + RV(29) * MDP(13) + RV(32) * MDP(14) + RV(36) * MDP(16) + RV(41) * MDP(17); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(14) * MDP(8) + RV(17) * MDP(9) + RV(20) * MDP(10) + RV(24) * MDP(11) + RV(27) * MDP(12) + RV(30) * MDP(13) + RV(33) * MDP(14) + RV(37) * MDP(16) + RV(42) * MDP(17); RV(12) * MDP(7) + RV(21) * MDP(10) + RV(38) * MDP(16) + RV(43) * MDP(17); RV(15) * MDP(8) + RV(18) * MDP(9) + RV(22) * MDP(10) + RV(39) * MDP(16) + RV(44) * MDP(17); RV(25) * MDP(11) + RV(28) * MDP(12) + RV(31) * MDP(13) + RV(34) * MDP(14) + RV(35) * MDP(15) + RV(40) * MDP(16) + RV(45) * MDP(17);];
tauJ = t1;
