% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:19
% EndTime: 2019-12-31 16:59:21
% DurationCPUTime: 0.52s
% Computational Cost: add. (306->129), mult. (752->170), div. (0->0), fcn. (303->2), ass. (0->67)
t114 = sin(qJ(2));
t134 = qJD(1) * t114;
t107 = pkin(5) * t134;
t93 = qJ(4) * t134 - t107;
t139 = qJD(3) - t93;
t146 = pkin(2) + pkin(3);
t85 = -t146 * qJD(2) + t139;
t111 = qJD(2) * qJ(3);
t115 = cos(qJ(2));
t133 = qJD(1) * t115;
t108 = pkin(5) * t133;
t95 = -qJ(4) * t133 + t108;
t90 = t111 + t95;
t112 = t114 ^ 2;
t113 = t115 ^ 2;
t149 = (t112 - t113) * MDP(5);
t125 = qJ(3) * t114 + pkin(1);
t98 = -pkin(2) * t115 - t125;
t91 = qJD(1) * t98;
t148 = t91 * MDP(13);
t147 = -0.2e1 * pkin(1);
t145 = pkin(5) - qJ(4);
t144 = pkin(5) * MDP(14);
t143 = qJD(2) * pkin(2);
t142 = qJ(3) * t115;
t92 = t146 * t115 + t125;
t141 = qJD(1) * t92;
t140 = t91 * MDP(11);
t82 = qJD(4) + t141;
t138 = qJD(4) + t82;
t128 = qJD(1) * qJD(2);
t126 = t115 * t128;
t131 = qJD(3) * t114;
t137 = qJ(3) * t126 + qJD(1) * t131;
t135 = qJ(4) * qJD(2);
t132 = qJD(2) * t114;
t130 = qJD(4) * t114;
t129 = qJD(4) * t115;
t101 = t145 * t115;
t127 = t114 * t128;
t79 = -t146 * t127 + t137;
t119 = -t146 * t114 + t142;
t81 = t119 * qJD(2) + t131;
t124 = qJD(1) * t81 + t79;
t84 = pkin(2) * t127 - t137;
t120 = pkin(2) * t114 - t142;
t89 = t120 * qJD(2) - t131;
t123 = -qJD(1) * t89 - t84;
t122 = MDP(12) + t144;
t121 = qJD(3) - t143;
t118 = qJD(1) ^ 2;
t117 = qJD(2) ^ 2;
t110 = qJD(2) * qJD(3);
t109 = 0.2e1 * t110;
t105 = pkin(5) * t126;
t103 = qJ(4) * t127;
t100 = t145 * t114;
t99 = t108 + t111;
t97 = t107 + t121;
t96 = -pkin(5) * t127 + t110;
t94 = t120 * qJD(1);
t88 = qJD(2) * t101 - t130;
t87 = -t145 * t132 - t129;
t86 = t119 * qJD(1);
t83 = t105 + (-t115 * t135 - t130) * qJD(1);
t80 = t103 + t110 + (-pkin(5) * t132 - t129) * qJD(1);
t1 = [(t84 * t98 + t89 * t91) * MDP(14) + (t100 * t83 + t101 * t80 + t79 * t92 + t81 * t82 + t85 * t88 + t87 * t90) * MDP(18) + (t123 * MDP(13) + t124 * MDP(16) + (-qJD(1) * t88 - t83) * MDP(17) + (-MDP(7) + (MDP(10) - MDP(13)) * pkin(5)) * t117) * t114 + (t117 * MDP(6) + t123 * MDP(11) + t96 * MDP(12) + t124 * MDP(15) + (-qJD(1) * t87 - t80) * MDP(17) + (t96 * MDP(14) + (-MDP(11) - MDP(9)) * t117) * pkin(5)) * t115 + (-t88 * MDP(15) + t87 * MDP(16) - 0.2e1 * qJD(1) * t149 + (qJD(1) * MDP(10) * t147 - 0.2e1 * t148 + (t82 + t141) * MDP(16) + (-qJD(1) * t100 - t85) * MDP(17) + t122 * t97) * t115 + (t140 - t82 * MDP(15) + t90 * MDP(17) - t122 * t99 + (t98 * MDP(11) - t92 * MDP(15) + t101 * MDP(17) + MDP(9) * t147 + (t122 * pkin(5) + (2 * MDP(4))) * t115) * qJD(1)) * t114) * qJD(2); t109 * MDP(13) + (qJ(3) * t96 + qJD(3) * t99 - t91 * t94) * MDP(14) + (qJD(2) * t95 - t105) * MDP(15) + (-qJD(2) * t93 + t103 + t109) * MDP(16) + (qJ(3) * t80 + t139 * t90 - t146 * t83 - t82 * t86 - t85 * t95) * MDP(18) + (-t114 * t115 * MDP(4) + t149 + (MDP(10) * t115 + MDP(9) * t114) * pkin(1)) * t118 + ((-t140 + (t99 - t111) * MDP(12) + t94 * MDP(13) + t99 * t144 + t138 * MDP(15) + (-pkin(5) * qJD(2) - t86) * MDP(16)) * t114 + (t94 * MDP(11) + (t121 - t97) * MDP(12) + t148 + (-t86 + t135) * MDP(15) - t138 * MDP(16) + (-t97 - t143) * t144) * t115) * qJD(1); (-qJD(2) * t99 + t105) * MDP(14) + (-qJ(4) * t126 - qJD(2) * t90 + t105) * MDP(18) + (MDP(13) + MDP(16)) * (-t112 * t118 - t117) + ((-MDP(11) - MDP(15)) * t118 * t115 + (t91 * MDP(14) - t138 * MDP(18)) * qJD(1)) * t114; t137 * MDP(18) + (-t112 - t113) * MDP(17) * t118 + ((t114 * t85 + t115 * t90) * MDP(18) + (0.2e1 * t115 * MDP(16) + (-t146 * MDP(18) - 0.2e1 * MDP(15)) * t114) * qJD(2)) * qJD(1);];
tauc = t1;
