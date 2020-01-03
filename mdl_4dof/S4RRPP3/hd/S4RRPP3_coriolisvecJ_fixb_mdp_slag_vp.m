% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:52
% EndTime: 2019-12-31 16:57:54
% DurationCPUTime: 0.55s
% Computational Cost: add. (569->139), mult. (1514->192), div. (0->0), fcn. (915->4), ass. (0->73)
t166 = (qJD(1) * qJD(2));
t181 = -2 * t166;
t148 = sin(qJ(2));
t180 = t148 * MDP(4);
t149 = cos(qJ(2));
t179 = (t148 ^ 2 - t149 ^ 2) * MDP(5);
t146 = sin(pkin(6));
t147 = cos(pkin(6));
t131 = t146 * t149 + t147 * t148;
t124 = t131 * qJD(1);
t120 = t124 ^ 2;
t177 = -qJ(3) - pkin(5);
t136 = t177 * t149;
t161 = t177 * t148;
t108 = -t136 * t146 - t147 * t161;
t158 = qJD(2) * t177;
t118 = qJD(3) * t149 + t148 * t158;
t113 = t118 * qJD(1);
t154 = -qJD(3) * t148 + t149 * t158;
t152 = t154 * qJD(1);
t94 = t113 * t146 - t147 * t152;
t176 = t94 * t108;
t134 = qJD(1) * t136;
t175 = t134 * t146;
t127 = t147 * t134;
t174 = t147 * t149;
t150 = qJD(2) ^ 2;
t173 = t148 * t150;
t172 = t149 * t150;
t95 = t147 * t113 + t146 * t152;
t133 = qJD(1) * t161;
t129 = qJD(2) * pkin(2) + t133;
t105 = t146 * t129 - t127;
t170 = qJD(1) * t149;
t169 = qJD(2) * t148;
t168 = t148 * qJD(1);
t107 = t133 * t147 + t175;
t167 = qJD(4) - t107;
t165 = qJD(2) * qJD(4);
t164 = pkin(2) * t169;
t163 = -pkin(2) * t149 - pkin(1);
t162 = t147 * t170;
t160 = t148 * t166;
t159 = t149 * t166;
t157 = pkin(1) * t181;
t156 = t163 * qJD(1);
t104 = t129 * t147 + t175;
t123 = t131 * qJD(2);
t115 = qJD(1) * t123;
t137 = t146 * t160;
t116 = t147 * t159 - t137;
t140 = pkin(2) * t160;
t155 = pkin(3) * t115 - qJ(4) * t116 + t140;
t135 = qJD(3) + t156;
t100 = t147 * t118 + t146 * t154;
t109 = -t147 * t136 + t146 * t161;
t121 = t146 * t168 - t162;
t99 = t118 * t146 - t147 * t154;
t153 = -t100 * t121 + t108 * t116 - t109 * t115 + t124 * t99 + t131 * t94;
t143 = -pkin(2) * t147 - pkin(3);
t141 = pkin(2) * t146 + qJ(4);
t130 = t146 * t148 - t174;
t126 = qJD(2) * t174 - t146 * t169;
t106 = t133 * t146 - t127;
t103 = pkin(3) * t130 - qJ(4) * t131 + t163;
t102 = qJD(2) * qJ(4) + t105;
t101 = -qJD(2) * pkin(3) + qJD(4) - t104;
t98 = pkin(2) * t168 + pkin(3) * t124 + qJ(4) * t121;
t97 = pkin(3) * t121 - qJ(4) * t124 + t135;
t93 = t165 + t95;
t92 = pkin(3) * t123 - qJ(4) * t126 - qJD(4) * t131 + t164;
t90 = -qJD(4) * t124 + t155;
t1 = [0.2e1 * t159 * t180 + t179 * t181 + MDP(6) * t172 - MDP(7) * t173 + (-pkin(5) * t172 + t148 * t157) * MDP(9) + (pkin(5) * t173 + t149 * t157) * MDP(10) + (-t104 * t126 - t105 * t123 - t130 * t95 + t153) * MDP(11) + (t105 * t100 - t104 * t99 + t176 + t95 * t109 + (t135 + t156) * t164) * MDP(12) + (-qJD(2) * t99 + t103 * t115 + t121 * t92 + t123 * t97 + t130 * t90) * MDP(13) + (t101 * t126 - t102 * t123 - t130 * t93 + t153) * MDP(14) + (qJD(2) * t100 - t103 * t116 - t124 * t92 - t97 * t126 - t90 * t131) * MDP(15) + (t100 * t102 + t101 * t99 + t103 * t90 + t109 * t93 + t92 * t97 + t176) * MDP(16); (t104 * t106 - t105 * t107) * MDP(12) + (qJD(2) * t106 - t94) * MDP(13) + (-t115 * t141 + t116 * t143) * MDP(14) + (-qJD(2) * t107 + 0.2e1 * t165 + t95) * MDP(15) + (-t101 * t106 + t167 * t102 + t141 * t93 + t143 * t94 - t97 * t98) * MDP(16) + ((t105 - t106) * MDP(11) - t97 * MDP(13) + (t102 - t106) * MDP(14) + t98 * MDP(15)) * t124 + ((-t104 + t107) * MDP(11) - t98 * MDP(13) + (t101 - t167) * MDP(14) - t97 * MDP(15)) * t121 + ((-t115 * t146 - t116 * t147) * MDP(11) + (-t135 * t168 + t146 * t95 - t147 * t94) * MDP(12)) * pkin(2) + (-t149 * t180 + t179 + (t149 * MDP(10) + t148 * MDP(9)) * pkin(1)) * qJD(1) ^ 2; (t104 * t124 + t105 * t121 + t140) * MDP(12) + t137 * MDP(15) + (t102 * t121 + (-qJD(4) - t101) * t124 + t155) * MDP(16) + ((t146 * t170 + t147 * t168 + t124) * MDP(13) + (t121 - t162) * MDP(15)) * qJD(2) + (MDP(11) + MDP(14)) * (-t121 ^ 2 - t120); t124 * t121 * MDP(13) + (-t137 + (t121 + t162) * qJD(2)) * MDP(14) + (-t120 - t150) * MDP(15) + (-qJD(2) * t102 + t124 * t97 + t94) * MDP(16);];
tauc = t1;
