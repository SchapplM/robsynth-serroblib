% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:56
% EndTime: 2019-12-31 17:36:57
% DurationCPUTime: 0.36s
% Computational Cost: add. (184->74), mult. (491->121), div. (0->0), fcn. (303->4), ass. (0->43)
t97 = sin(pkin(8));
t95 = t97 ^ 2;
t98 = cos(pkin(8));
t128 = t98 ^ 2 + t95;
t107 = t97 * qJ(4) + pkin(2);
t127 = qJD(2) * (t98 * pkin(3) + t107);
t126 = t128 * (MDP(7) + MDP(10));
t114 = qJ(3) * qJD(2);
t124 = (t97 * qJD(1) + t114 * t98) * t98;
t99 = sin(qJ(5));
t123 = t98 * t99;
t122 = -pkin(6) + qJ(3);
t100 = cos(qJ(5));
t120 = t97 * t100;
t119 = t98 * MDP(9);
t80 = t98 * t100 + t97 * t99;
t118 = qJD(2) * t80;
t117 = qJD(3) * t97;
t116 = qJD(4) * t97;
t115 = t95 * MDP(11);
t113 = qJD(2) * qJD(3);
t112 = t128 * qJ(3) * t113 + qJD(3) * t124;
t111 = qJD(2) * t123;
t109 = qJD(2) * t120;
t108 = qJD(5) * t120;
t82 = t98 * qJD(1) - t97 * t114;
t81 = t120 - t123;
t104 = t80 * qJD(3);
t72 = t80 * qJD(5);
t103 = t81 * qJD(3);
t77 = (pkin(3) + pkin(4)) * t98 + t107;
t102 = qJD(2) ^ 2;
t87 = qJD(5) * t111;
t86 = t122 * t98;
t85 = t122 * t97;
t79 = qJD(4) - t82;
t76 = t109 - t111;
t73 = -qJD(5) * t123 + t108;
t70 = qJD(3) - t127;
t68 = qJD(2) * t108 - t87;
t67 = qJD(2) * t72;
t65 = qJD(2) * t77 - qJD(3);
t1 = [(-t73 * MDP(18) + t72 * MDP(19)) * qJD(5); (-t117 * t82 + t112) * MDP(8) + (t79 * t117 + (-t70 + t127) * t116 + t112) * MDP(12) + (-t67 * t81 - t76 * t72) * MDP(13) + (t118 * t72 + t67 * t80 - t81 * t68 - t76 * t73) * MDP(14) - t72 * qJD(5) * MDP(15) - t73 * qJD(5) * MDP(16) + (t65 * t73 + t77 * t68 + ((-t100 * t86 - t85 * t99) * qJD(5) + t103) * qJD(5) + 0.2e1 * t118 * t116) * MDP(18) + (-t65 * t72 - t77 * t67 + ((-t100 * t85 + t86 * t99) * qJD(5) - t104) * qJD(5) + (qJD(2) * t81 + t76) * t116) * MDP(19) + 0.2e1 * t113 * t126 + 0.2e1 * (t97 * t119 + t115) * qJD(2) * qJD(4); (-t76 * qJD(5) + t87) * MDP(18) + t118 * qJD(5) * MDP(19) - t102 * t126 + ((t82 * t97 - t124) * MDP(8) + (-t79 * t97 - t116 - t124) * MDP(12) + (-MDP(18) * t120 + MDP(19) * t80) * qJD(5)) * qJD(2); -t102 * t115 + (-t99 * MDP(18) - t100 * MDP(19)) * qJD(5) ^ 2 + (-t102 * t119 + ((qJD(3) + t70) * MDP(12) - t118 * MDP(18) - t76 * MDP(19)) * qJD(2)) * t97; t76 * t118 * MDP(13) + (-t118 ^ 2 + t76 ^ 2) * MDP(14) + (t87 + (t76 - t109) * qJD(5)) * MDP(16) + (qJD(2) * t103 - t65 * t76) * MDP(18) + (-qJD(2) * t104 + t118 * t65) * MDP(19);];
tauc = t1;
