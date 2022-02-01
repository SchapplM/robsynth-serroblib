% Calculate vector of inverse dynamics joint torques for
% S5RPPRP1
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [49x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5RPPRP1_invdynJ_fixb_regmin2vec.m
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRP1_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(49,1), zeros(18,1)}
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vr: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:51
% EndTime: 2022-01-23 09:12:51
% DurationCPUTime: 0.03s
% Computational Cost: add. (44->44), mult. (49->49), div. (0->0), fcn. (49->49), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(8) + RV(15) * MDP(9) + RV(17) * MDP(10) + RV(19) * MDP(11) + RV(21) * MDP(12) + RV(23) * MDP(13) + RV(27) * MDP(14) + RV(31) * MDP(15) + RV(36) * MDP(16) + RV(41) * MDP(17) + RV(45) * MDP(18); RV(5) * MDP(4) + RV(11) * MDP(7) + RV(24) * MDP(13) + RV(28) * MDP(14) + RV(32) * MDP(15) + RV(37) * MDP(16) + RV(46) * MDP(18); RV(7) * MDP(5) + RV(9) * MDP(6) + RV(12) * MDP(7) + RV(25) * MDP(13) + RV(29) * MDP(14) + RV(33) * MDP(15) + RV(38) * MDP(16) + RV(42) * MDP(17) + RV(47) * MDP(18); RV(14) * MDP(8) + RV(16) * MDP(9) + RV(18) * MDP(10) + RV(20) * MDP(11) + RV(22) * MDP(12) + RV(26) * MDP(13) + RV(30) * MDP(14) + RV(34) * MDP(15) + RV(39) * MDP(16) + RV(43) * MDP(17) + RV(48) * MDP(18); RV(35) * MDP(15) + RV(40) * MDP(16) + RV(44) * MDP(17) + RV(49) * MDP(18);];
tauJ = t1;
