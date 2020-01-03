% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPP5
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
%   see S4RRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:38
% EndTime: 2019-12-31 17:00:40
% DurationCPUTime: 0.54s
% Computational Cost: add. (306->125), mult. (743->160), div. (0->0), fcn. (294->2), ass. (0->70)
t120 = cos(qJ(2));
t139 = qJD(1) * t120;
t110 = pkin(5) * t139;
t96 = pkin(3) * t139 + t110;
t156 = -qJD(4) - t96;
t119 = sin(qJ(2));
t150 = pkin(3) + pkin(5);
t153 = t150 * t119;
t155 = qJD(1) * t153;
t138 = qJD(2) * qJ(3);
t90 = t138 - t156;
t116 = t119 ^ 2;
t117 = t120 ^ 2;
t154 = (t116 - t117) * MDP(5);
t118 = -pkin(2) - qJ(4);
t152 = qJD(2) * t118;
t149 = MDP(14) * pkin(5);
t127 = MDP(11) + t149;
t151 = -2 * pkin(1);
t148 = qJD(2) * pkin(2);
t131 = -t119 * qJ(3) - pkin(1);
t92 = t118 * t120 + t131;
t82 = qJD(1) * t92;
t95 = qJD(2) * t153;
t147 = qJD(1) * t95;
t122 = qJD(1) ^ 2;
t146 = t120 * t122;
t145 = t82 * MDP(17);
t100 = -pkin(2) * t120 + t131;
t91 = qJD(1) * t100;
t144 = t91 * MDP(13);
t140 = qJD(1) * t119;
t109 = pkin(5) * t140;
t94 = -pkin(3) * t140 - t109;
t143 = qJD(3) - t94;
t136 = qJD(1) * qJD(2);
t132 = t120 * t136;
t105 = pkin(5) * t132;
t142 = pkin(3) * t132 + t105;
t137 = t119 * qJD(3);
t135 = qJD(2) * qJD(4);
t103 = t150 * t120;
t133 = t119 * t136;
t106 = pkin(2) * t133;
t125 = -qJ(3) * t120 + qJ(4) * t119;
t123 = t125 * qJD(2) - t120 * qJD(4) - t137;
t79 = t123 * qJD(1) + t106;
t112 = t119 * t148;
t80 = t112 + t123;
t129 = -qJD(1) * t80 - t79;
t124 = -t120 * t138 - t137;
t83 = t124 * qJD(1) + t106;
t89 = t112 + t124;
t128 = qJD(1) * t89 + t83;
t87 = -t135 + t142;
t121 = qJD(2) ^ 2;
t115 = qJD(2) * qJD(3);
t114 = 0.2e1 * t115;
t113 = pkin(2) * t140;
t101 = -t110 - t138;
t99 = qJD(3) + t109 - t148;
t98 = pkin(5) * t133 - t115;
t97 = qJD(2) * t103;
t93 = -qJ(3) * t139 + t113;
t88 = t125 * qJD(1) + t113;
t86 = t115 - t147;
t85 = t143 + t152;
t84 = t91 * t140;
t81 = t82 * t139;
t1 = [(t100 * t83 + t89 * t91) * MDP(14) + (t103 * t86 + t153 * t87 + t79 * t92 + t80 * t82 + t85 * t97 - t90 * t95) * MDP(18) + (-t128 * MDP(13) + (qJD(1) * t97 + t87) * MDP(15) + t129 * MDP(16) + (-MDP(7) + (MDP(10) - MDP(13)) * pkin(5)) * t121) * t119 + (t121 * MDP(6) - t98 * MDP(11) + t128 * MDP(12) + (t86 - t147) * MDP(15) + t129 * MDP(17) + (-t98 * MDP(14) + (MDP(12) - MDP(9)) * t121) * pkin(5)) * t120 + (-t95 * MDP(16) - t97 * MDP(17) - 0.2e1 * qJD(1) * t154 + (qJD(1) * MDP(10) * t151 - 0.2e1 * t144 + (t85 + t155) * MDP(15) - 0.2e1 * t82 * MDP(16) + t127 * t99) * t120 + (-t91 * MDP(12) - t90 * MDP(15) + t145 + t127 * t101 + (-t100 * MDP(12) - t103 * MDP(15) + t92 * MDP(17) + MDP(9) * t151 + (t127 * pkin(5) + (2 * MDP(4))) * t120) * qJD(1)) * t119) * qJD(2); t84 * MDP(12) + t114 * MDP(13) + (-qJ(3) * t98 - qJD(3) * t101 - t91 * t93) * MDP(14) + (-qJD(2) * t94 + t114 + t81) * MDP(16) + (qJD(2) * t96 + 0.2e1 * t135 - t142) * MDP(17) + (qJ(3) * t86 + t118 * t87 + t143 * t90 + t156 * t85 - t82 * t88) * MDP(18) + (MDP(11) + MDP(15)) * qJD(3) * t139 + (-t119 * t120 * MDP(4) + t154 + (MDP(10) * t120 + MDP(9) * t119) * pkin(1)) * t122 + ((-t93 * MDP(12) + t144 + (-t85 - t94 + t152) * MDP(15) + t88 * MDP(17) + t127 * (-t99 - t148)) * t120 + ((-t101 - t138) * MDP(11) + t93 * MDP(13) - t101 * t149 + (-t150 * qJD(2) + t88) * MDP(16) - t145) * t119) * qJD(1); (qJD(2) * t101 + t105 + t84) * MDP(14) + (-qJD(2) * t90 + t87) * MDP(18) + (qJD(1) * t82 * MDP(18) + (MDP(12) - MDP(17)) * t146) * t119 + (MDP(13) + MDP(16)) * (-t116 * t122 - t121); -t119 * MDP(16) * t146 + (-t117 * t122 - t121) * MDP(17) + (t115 + t81 + (t85 - t155) * qJD(2)) * MDP(18);];
tauc = t1;
