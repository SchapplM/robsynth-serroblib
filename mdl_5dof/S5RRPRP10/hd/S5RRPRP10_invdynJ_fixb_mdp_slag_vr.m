% Calculate vector of inverse dynamics joint torques for
% S5RRPRP10
% The function exploits the sparsity of the regressor matrix
% 
% Input:
% RV [71x1]
%   vector of non-Null entries of the regressor matrix. (columns, then rows).
%   see S5RRPRP10_invdynJ_fixb_regmin2vec.m
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP10_invdynJ_fixb_mdp_slag_vr(RV, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(71,1), zeros(25,1)}
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vr: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_mult_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:03:07
% EndTime: 2021-01-15 21:03:07
% DurationCPUTime: 0.03s
% Computational Cost: add. (66->66), mult. (71->71), div. (0->0), fcn. (71->71), ass. (0->1)
t1 = [RV(1) * MDP(1) + RV(2) * MDP(2) + RV(3) * MDP(3) + RV(4) * MDP(4) + RV(6) * MDP(5) + RV(8) * MDP(6) + RV(10) * MDP(7) + RV(13) * MDP(9) + RV(15) * MDP(10) + RV(17) * MDP(11) + RV(20) * MDP(12) + RV(23) * MDP(13) + RV(26) * MDP(14) + RV(29) * MDP(15) + RV(32) * MDP(16) + RV(35) * MDP(17) + RV(38) * MDP(18) + RV(41) * MDP(19) + RV(44) * MDP(20) + RV(48) * MDP(21) + RV(52) * MDP(22) + RV(57) * MDP(23) + RV(62) * MDP(24) + RV(67) * MDP(25); RV(5) * MDP(4) + RV(7) * MDP(5) + RV(9) * MDP(6) + RV(11) * MDP(7) + RV(12) * MDP(8) + RV(14) * MDP(9) + RV(16) * MDP(10) + RV(18) * MDP(11) + RV(21) * MDP(12) + RV(24) * MDP(13) + RV(27) * MDP(14) + RV(30) * MDP(15) + RV(33) * MDP(16) + RV(36) * MDP(17) + RV(39) * MDP(18) + RV(42) * MDP(19) + RV(45) * MDP(20) + RV(49) * MDP(21) + RV(53) * MDP(22) + RV(58) * MDP(23) + RV(63) * MDP(24) + RV(68) * MDP(25); RV(19) * MDP(11) + RV(22) * MDP(12) + RV(25) * MDP(13) + RV(28) * MDP(14) + RV(46) * MDP(20) + RV(50) * MDP(21) + RV(54) * MDP(22) + RV(59) * MDP(23) + RV(64) * MDP(24) + RV(69) * MDP(25); RV(31) * MDP(15) + RV(34) * MDP(16) + RV(37) * MDP(17) + RV(40) * MDP(18) + RV(43) * MDP(19) + RV(47) * MDP(20) + RV(51) * MDP(21) + RV(55) * MDP(22) + RV(60) * MDP(23) + RV(65) * MDP(24) + RV(70) * MDP(25); RV(56) * MDP(22) + RV(61) * MDP(23) + RV(66) * MDP(24) + RV(71) * MDP(25);];
tauJ = t1;
