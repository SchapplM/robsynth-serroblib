% Calculate vector of inverse dynamics joint torques for
% S5RPRPR14
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [65x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5RPRPR14_invdynJ_fixb_regmin2vec.m
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR14_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(65,1), zeros(24,1)}
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vr: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:21
% EndTime: 2021-01-15 12:17:21
% DurationCPUTime: 0.03s
% Computational Cost: add. (60->60), mult. (65->65), div. (0->0), fcn. (65->65), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(12) * MDP(8) + RV(14) * MDP(9) + RV(16) * MDP(10) + RV(19) * MDP(12) + RV(22) * MDP(13) + RV(25) * MDP(14) + RV(29) * MDP(15) + RV(33) * MDP(16) + RV(37) * MDP(17) + RV(41) * MDP(18) + RV(44) * MDP(19) + RV(47) * MDP(20) + RV(50) * MDP(21) + RV(53) * MDP(22) + RV(56) * MDP(23) + RV(61) * MDP(24); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(20) * MDP(12) + RV(23) * MDP(13) + RV(26) * MDP(14) + RV(30) * MDP(15) + RV(34) * MDP(16) + RV(38) * MDP(17) + RV(57) * MDP(23) + RV(62) * MDP(24); RV(11) * MDP(7) + RV(13) * MDP(8) + RV(15) * MDP(9) + RV(17) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(12) + RV(24) * MDP(13) + RV(27) * MDP(14) + RV(31) * MDP(15) + RV(35) * MDP(16) + RV(39) * MDP(17) + RV(42) * MDP(18) + RV(45) * MDP(19) + RV(48) * MDP(20) + RV(51) * MDP(21) + RV(54) * MDP(22) + RV(58) * MDP(23) + RV(63) * MDP(24); RV(28) * MDP(14) + RV(32) * MDP(15) + RV(36) * MDP(16) + RV(40) * MDP(17) + RV(59) * MDP(23) + RV(64) * MDP(24); RV(43) * MDP(18) + RV(46) * MDP(19) + RV(49) * MDP(20) + RV(52) * MDP(21) + RV(55) * MDP(22) + RV(60) * MDP(23) + RV(65) * MDP(24);];
tauJ = t1;
