% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:15
% EndTime: 2019-12-05 16:40:19
% DurationCPUTime: 0.57s
% Computational Cost: add. (445->109), mult. (800->166), div. (0->0), fcn. (370->4), ass. (0->68)
t164 = qJ(5) + pkin(7);
t116 = sin(qJ(4));
t163 = MDP(8) * t116;
t118 = cos(qJ(4));
t113 = qJD(2) + qJD(3);
t117 = sin(qJ(3));
t156 = pkin(2) * qJD(2);
t139 = t117 * t156;
t129 = t164 * t113 + t139;
t88 = t116 * qJD(1) + t129 * t118;
t162 = qJD(4) * t88;
t114 = t116 ^ 2;
t115 = t118 ^ 2;
t161 = (t114 - t115) * MDP(9);
t87 = t118 * qJD(1) - t129 * t116;
t160 = MDP(13) * t116 + MDP(14) * t118;
t119 = cos(qJ(3));
t159 = pkin(2) * t119;
t154 = qJD(4) * pkin(4);
t84 = t87 + t154;
t157 = t84 - t87;
t155 = pkin(2) * qJD(3);
t136 = t119 * t156;
t100 = -pkin(3) * t113 - t136;
t135 = qJD(2) * t155;
t105 = t117 * t135;
t143 = qJD(4) * t118;
t153 = t100 * t143 + t116 * t105;
t152 = t113 * t116;
t120 = qJD(4) ^ 2;
t151 = t116 * t120;
t150 = t118 * t120;
t107 = pkin(2) * t117 + pkin(7);
t149 = -qJ(5) - t107;
t144 = qJD(4) * t116;
t133 = t113 * t144;
t93 = pkin(4) * t133 + t105;
t148 = -t114 - t115;
t142 = qJD(4) * t119;
t140 = -qJD(3) + t113;
t138 = t119 * t155;
t137 = pkin(4) * t144;
t134 = -pkin(4) * t118 - pkin(3);
t132 = t113 * t143;
t127 = t119 * t135;
t121 = qJD(5) * t113 + t127;
t82 = t87 * qJD(4) + t121 * t118;
t83 = -t121 * t116 - t162;
t131 = -t83 * t116 + t82 * t118;
t130 = qJD(4) * t164;
t128 = qJD(4) * t149;
t126 = qJD(4) * t138;
t125 = t116 * t84 - t118 * t88;
t123 = -0.2e1 * t113 * qJD(4) * t161 + MDP(10) * t150 - MDP(11) * t151 + 0.2e1 * t132 * t163;
t112 = t113 ^ 2;
t111 = t118 * qJ(5);
t109 = t118 * qJD(5);
t103 = pkin(7) * t118 + t111;
t102 = t164 * t116;
t97 = t107 * t118 + t111;
t96 = t149 * t116;
t94 = t100 * t144;
t92 = -t116 * qJD(5) - t118 * t130;
t91 = -t116 * t130 + t109;
t90 = t134 * t113 + qJD(5) - t136;
t86 = (-qJD(5) - t138) * t116 + t118 * t128;
t85 = t116 * t128 + t118 * t138 + t109;
t1 = [(-t120 * MDP(14) + (t83 + t162) * MDP(16)) * t118 + (-t120 * MDP(13) + (-qJD(4) * t84 + t82) * MDP(16)) * t116; -t105 * MDP(6) - MDP(7) * t127 + (-t118 * t105 - t107 * t150 - t116 * t126 + t94) * MDP(13) + (t107 * t151 - t118 * t126 + t153) * MDP(14) + (-t84 * t143 - t88 * t144 + t131) * MDP(15) + (t82 * t97 + t88 * t85 + t83 * t96 + t84 * t86 + t93 * (t134 - t159) + t90 * (t117 * t155 + t137)) * MDP(16) + ((-t116 * t86 + t118 * t85) * MDP(15) + ((-t116 * t97 - t118 * t96) * MDP(15) + t160 * (-pkin(3) - t159)) * qJD(4) + (-t119 * MDP(7) + (-MDP(13) * t118 + MDP(14) * t116 - MDP(6)) * t117) * t155) * t113 + t123; (t113 * t139 - t105) * MDP(6) + t140 * MDP(7) * t136 + (-pkin(3) * t133 - pkin(7) * t150 + t94 + (t140 * t118 * t117 + t116 * t142) * t156) * MDP(13) + (-pkin(3) * t132 + pkin(7) * t151 + (-t117 * t152 + t118 * t142) * t156 + t153) * MDP(14) + ((-t116 * t88 - t118 * t84) * qJD(4) + (-t116 * t92 + t118 * t91 + (t102 * t118 - t103 * t116) * qJD(4) + t148 * t136) * t113 + t131) * MDP(15) + (t82 * t103 + t88 * t91 - t83 * t102 + t84 * t92 + t93 * t134 + t90 * t137 + (-t117 * t90 + t125 * t119) * t156) * MDP(16) + t123; (-t154 + t157) * t118 * t113 * MDP(15) + (t157 * t88 + (-t90 * t152 + t83) * pkin(4)) * MDP(16) + t160 * (-t100 * t113 - t127) + (-t118 * t163 + t161) * t112; (t125 * t113 + t93) * MDP(16) + t148 * MDP(15) * t112;];
tauc = t1;
