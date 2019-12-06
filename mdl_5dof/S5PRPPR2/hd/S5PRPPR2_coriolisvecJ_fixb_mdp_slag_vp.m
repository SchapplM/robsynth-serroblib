% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:43
% EndTime: 2019-12-05 15:24:46
% DurationCPUTime: 0.56s
% Computational Cost: add. (343->96), mult. (908->153), div. (0->0), fcn. (671->8), ass. (0->58)
t116 = sin(pkin(9));
t118 = cos(pkin(9));
t138 = t116 ^ 2 + t118 ^ 2;
t117 = sin(pkin(8));
t121 = sin(qJ(2));
t137 = qJD(1) * t121;
t108 = t117 * t137;
t119 = cos(pkin(8));
t123 = cos(qJ(2));
t136 = qJD(1) * t123;
t131 = t119 * t136;
t93 = -t108 + t131;
t144 = qJD(4) - t93;
t122 = cos(qJ(5));
t120 = sin(qJ(5));
t139 = t116 * t120;
t100 = -t118 * t122 + t139;
t95 = t100 * qJD(5);
t102 = t116 * t122 + t118 * t120;
t94 = t102 * qJD(2);
t143 = t138 * qJD(2);
t142 = pkin(2) * t119;
t101 = t117 * t123 + t119 * t121;
t88 = t101 * qJD(2);
t84 = qJD(1) * t88;
t99 = t117 * t121 - t119 * t123;
t141 = t84 * t99;
t111 = pkin(2) * t117 + qJ(4);
t140 = pkin(6) + t111;
t107 = qJD(2) * pkin(2) + t136;
t109 = t119 * t137;
t83 = t117 * t107 + t109;
t135 = qJD(2) * t116;
t134 = qJD(2) * t118;
t132 = -pkin(4) * t118 - pkin(3);
t130 = t120 * t135;
t129 = t122 * t134;
t105 = qJD(2) * t131;
t81 = t105 + (qJD(4) - t108) * qJD(2);
t128 = t138 * t81;
t82 = t107 * t119 - t108;
t126 = qJD(4) - t82;
t125 = t138 * (qJD(2) * qJ(4) + t83);
t96 = t102 * qJD(5);
t124 = qJD(2) ^ 2;
t106 = qJD(5) * t129;
t104 = t132 - t142;
t98 = t140 * t118;
t97 = t140 * t116;
t92 = t99 * qJD(2);
t90 = -t129 + t130;
t89 = t117 * t136 + t109;
t87 = qJD(2) * t96;
t86 = -qJD(5) * t130 + t106;
t85 = -qJD(2) * t108 + t105;
t79 = -qJD(2) * pkin(3) + t126;
t77 = t132 * qJD(2) + t126;
t1 = [(t101 * t85 - t83 * t92 + t141) * MDP(5) - t92 * MDP(8) * t143 + (t101 * t128 - t125 * t92 + t141) * MDP(9) + (t99 * t87 + (t101 * t95 + t102 * t92) * qJD(5)) * MDP(15) + (t99 * t86 + (-t100 * t92 + t101 * t96) * qJD(5)) * MDP(16) + (t90 * MDP(15) + t94 * MDP(16) - t82 * MDP(5) - MDP(6) * t134 + MDP(7) * t135 + t79 * MDP(9)) * t88 + (-t121 * MDP(3) - t123 * MDP(4)) * t124; (t82 * t89 - t83 * t93 + (t117 * t85 - t119 * t84) * pkin(2)) * MDP(5) + (t144 * t143 + t128) * MDP(8) + (t84 * (-pkin(3) - t142) - t79 * t89 + t111 * t128 + t144 * t125) * MDP(9) + (t102 * t86 - t94 * t95) * MDP(10) + (-t100 * t86 - t102 * t87 + t90 * t95 - t94 * t96) * MDP(11) - t95 * qJD(5) * MDP(12) - t96 * qJD(5) * MDP(13) + (t84 * t100 + t104 * t87 + t77 * t96 - t89 * t90 + ((t120 * t97 - t122 * t98) * qJD(5) - t144 * t102) * qJD(5)) * MDP(15) + (t84 * t102 + t104 * t86 - t77 * t95 - t89 * t94 + ((t120 * t98 + t122 * t97) * qJD(5) + t144 * t100) * qJD(5)) * MDP(16) + (t118 * MDP(6) - t116 * MDP(7)) * (qJD(2) * t89 - t84); (-t96 * MDP(15) + t95 * MDP(16)) * qJD(5); t94 * qJD(5) * MDP(15) + (-qJD(5) * t90 + t106) * MDP(16) - t138 * MDP(8) * t124 + ((t101 * qJD(1) - t125) * MDP(9) + (t102 * MDP(15) - MDP(16) * t139) * qJD(5)) * qJD(2); t94 * t90 * MDP(10) + (-t90 ^ 2 + t94 ^ 2) * MDP(11) + (t106 + (t90 - t130) * qJD(5)) * MDP(12) + (-t102 * t81 - t77 * t94) * MDP(15) + (t100 * t81 + t77 * t90) * MDP(16);];
tauc = t1;
