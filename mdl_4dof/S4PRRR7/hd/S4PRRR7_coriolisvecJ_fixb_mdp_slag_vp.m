% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:50
% EndTime: 2019-12-31 16:36:53
% DurationCPUTime: 1.10s
% Computational Cost: add. (428->168), mult. (1177->277), div. (0->0), fcn. (808->8), ass. (0->84)
t131 = sin(qJ(3));
t189 = MDP(5) * t131;
t126 = t131 ^ 2;
t134 = cos(qJ(3));
t188 = (-t134 ^ 2 + t126) * MDP(6);
t187 = t134 * MDP(11);
t186 = qJD(2) * pkin(2);
t132 = sin(qJ(2));
t128 = sin(pkin(4));
t170 = qJD(1) * t128;
t155 = t132 * t170;
t119 = qJD(2) * pkin(6) + t155;
t129 = cos(pkin(4));
t135 = cos(qJ(2));
t177 = t128 * t135;
t153 = qJD(2) * t177;
t140 = qJD(1) * (qJD(3) * t129 + t153);
t165 = qJD(3) * t134;
t100 = t119 * t165 + t131 * t140;
t130 = sin(qJ(4));
t185 = t100 * t130;
t133 = cos(qJ(4));
t184 = t100 * t133;
t164 = qJD(4) * t130;
t150 = t131 * t164;
t157 = qJD(2) * qJD(3);
t148 = t134 * t157;
t156 = qJD(3) * qJD(4);
t172 = (t148 + t156) * t133;
t106 = -qJD(2) * t150 + t172;
t183 = t106 * t130;
t160 = t130 * qJD(2);
t152 = t131 * t160;
t158 = t133 * qJD(3);
t114 = t152 - t158;
t168 = qJD(2) * t134;
t123 = -qJD(4) + t168;
t182 = t114 * t123;
t159 = t133 * qJD(2);
t167 = qJD(3) * t130;
t116 = t131 * t159 + t167;
t181 = t116 * t123;
t180 = t123 * t130;
t179 = t123 * t133;
t178 = t128 * t132;
t176 = t130 * t134;
t136 = qJD(3) ^ 2;
t175 = t131 * t136;
t174 = t133 * t134;
t173 = t134 * t136;
t169 = qJD(1) * t129;
t166 = qJD(3) * t131;
t163 = qJD(4) * t133;
t162 = qJD(4) * t135;
t109 = -t119 * t131 + t134 * t169;
t104 = -qJD(3) * pkin(3) - t109;
t161 = t104 * qJD(4);
t154 = t135 * t170;
t151 = t123 * t163;
t149 = t131 * t163;
t147 = MDP(16) * t166;
t146 = pkin(3) * t131 - pkin(7) * t134;
t118 = t146 * qJD(3);
t145 = -t118 + t155;
t110 = t119 * t134 + t131 * t169;
t105 = qJD(3) * pkin(7) + t110;
t121 = -pkin(3) * t134 - pkin(7) * t131 - pkin(2);
t111 = t121 * qJD(2) - t154;
t98 = t105 * t133 + t111 * t130;
t143 = t105 * t130 - t111 * t133;
t142 = qJD(2) * t126 - t123 * t134;
t113 = t129 * t131 + t134 * t178;
t112 = -t129 * t134 + t131 * t178;
t139 = -0.2e1 * qJD(3) * t186;
t99 = -t119 * t166 + t134 * t140;
t138 = qJD(3) * t104 + qJD(4) * t111 - t123 * t154 + t99;
t137 = qJD(2) ^ 2;
t117 = t146 * qJD(2);
t108 = (t118 + t155) * qJD(2);
t107 = t130 * t156 + (t130 * t165 + t149) * qJD(2);
t103 = t133 * t108;
t102 = t113 * qJD(3) + t131 * t153;
t101 = -t112 * qJD(3) + t134 * t153;
t1 = [(-(-t101 * t130 - t113 * t163) * t123 + t102 * t114 + t112 * t107) * MDP(17) + ((t101 * t133 - t113 * t164) * t123 + t102 * t116 + t112 * t106) * MDP(18) + ((-(t130 * t162 + t132 * t159) * MDP(17) + (t132 * t160 - t133 * t162) * MDP(18)) * t123 + (-t135 * MDP(4) + (-MDP(10) * t134 + MDP(11) * t131 - MDP(3)) * t132) * t137) * t128 + (-t102 * MDP(10) - t101 * MDP(11) + (-t177 * t187 + (-MDP(10) * t177 + (-t113 * t130 - t133 * t177) * MDP(17) - (t113 * t133 - t130 * t177) * MDP(18)) * t131) * qJD(2)) * qJD(3); 0.2e1 * t148 * t189 - 0.2e1 * t157 * t188 + MDP(7) * t173 - MDP(8) * t175 + (-pkin(6) * t173 + t131 * t139) * MDP(10) + (pkin(6) * t175 + t134 * t139) * MDP(11) + (t106 * t131 * t133 + (t134 * t158 - t150) * t116) * MDP(12) + ((-t114 * t133 - t116 * t130) * t165 + (-t183 - t107 * t133 + (t114 * t130 - t116 * t133) * qJD(4)) * t131) * MDP(13) + (t123 * t150 - t106 * t134 + (t116 * t131 + t142 * t133) * qJD(3)) * MDP(14) + (t123 * t149 + t107 * t134 + (-t114 * t131 - t142 * t130) * qJD(3)) * MDP(15) + (-t123 - t168) * t147 + ((t121 * t164 + t145 * t133) * t123 + (t105 * t163 - t103 + (qJD(3) * t114 + t151) * pkin(6) + t138 * t130) * t134 + (-t114 * t154 + t133 * t161 + pkin(6) * t107 + t185 + (-pkin(6) * t180 + (-pkin(6) * t176 + t121 * t133) * qJD(2) - t143) * qJD(3)) * t131) * MDP(17) + ((t121 * t163 - t145 * t130) * t123 + (qJD(3) * pkin(6) * t116 + (t108 + (-pkin(6) * t123 - t105) * qJD(4)) * t130 + t138 * t133) * t134 + (-t116 * t154 - t130 * t161 + pkin(6) * t106 + t184 + (-pkin(6) * t179 - (pkin(6) * t174 + t121 * t130) * qJD(2) - t98) * qJD(3)) * t131) * MDP(18); (-t116 * t179 + t183) * MDP(12) + ((t106 + t182) * t133 + (-t107 + t181) * t130) * MDP(13) - t151 * MDP(14) + t123 * t164 * MDP(15) + (-pkin(3) * t107 - t184 + (-t109 * t130 + t117 * t133) * t123 - t110 * t114 + (pkin(7) * t179 + t104 * t130) * qJD(4)) * MDP(17) + (-pkin(3) * t106 + t185 - (t109 * t133 + t117 * t130) * t123 - t110 * t116 + (-pkin(7) * t180 + t104 * t133) * qJD(4)) * MDP(18) + (-t134 * t189 + t188) * t137 + ((t123 * t174 + (-t116 + t167) * t131) * MDP(14) + (-t123 * t176 + (t114 + t158) * t131) * MDP(15) + t123 * t131 * MDP(16) + (t143 * t131 + (-pkin(7) * t166 - t104 * t134) * t130) * MDP(17) + (-t104 * t174 + (-pkin(7) * t158 + t98) * t131) * MDP(18) + (t131 * MDP(10) + t187) * t186) * qJD(2); t116 * t114 * MDP(12) + (-t114 ^ 2 + t116 ^ 2) * MDP(13) + (t172 - t182) * MDP(14) + (-t130 * t148 - t181) * MDP(15) + qJD(2) * t147 + (-t104 * t116 - t123 * t98 - t130 * t99 + t103) * MDP(17) + (t104 * t114 - t130 * t108 + t123 * t143 - t133 * t99) * MDP(18) + (-MDP(14) * t152 - t116 * MDP(15) - t98 * MDP(17) + t143 * MDP(18)) * qJD(4);];
tauc = t1;
