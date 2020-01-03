% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:58
% EndTime: 2019-12-31 17:12:00
% DurationCPUTime: 1.13s
% Computational Cost: add. (470->185), mult. (1162->276), div. (0->0), fcn. (610->4), ass. (0->97)
t185 = (qJD(1) * qJD(2));
t214 = -2 * t185;
t156 = sin(qJ(2));
t153 = t156 ^ 2;
t158 = cos(qJ(2));
t154 = t158 ^ 2;
t213 = (t153 - t154) * MDP(5);
t211 = pkin(3) + pkin(5);
t212 = t211 * t156;
t159 = -pkin(2) - pkin(6);
t210 = qJD(2) * pkin(2);
t157 = cos(qJ(4));
t155 = sin(qJ(4));
t193 = qJD(2) * t155;
t195 = qJD(1) * t158;
t126 = t157 * t195 + t193;
t179 = t156 * t185;
t140 = t155 * t179;
t111 = -t126 * qJD(4) + t140;
t209 = t111 * t157;
t118 = (-qJD(1) * t212 + qJD(3)) * qJD(2);
t208 = t118 * t155;
t207 = t118 * t157;
t196 = qJD(1) * t156;
t147 = qJD(4) + t196;
t206 = t126 * t147;
t180 = t155 * t195;
t191 = qJD(2) * t157;
t128 = -t180 + t191;
t205 = t128 * t147;
t204 = t128 * t158;
t203 = t155 * t147;
t160 = qJD(2) ^ 2;
t202 = t156 * t160;
t201 = t157 * t147;
t200 = t158 * t160;
t161 = qJD(1) ^ 2;
t199 = t158 * t161;
t198 = qJD(4) * t180 + t157 * t179;
t176 = -qJ(3) * t156 - pkin(1);
t136 = -pkin(2) * t158 + t176;
t123 = qJD(1) * t136;
t194 = qJD(2) * qJ(3);
t192 = qJD(2) * t156;
t190 = qJD(2) * t158;
t189 = qJD(3) * t156;
t188 = qJD(4) * t157;
t149 = pkin(5) * t195;
t137 = -t149 - t194;
t150 = pkin(3) * t195;
t122 = -t137 + t150;
t187 = t122 * qJD(4);
t148 = pkin(5) * t196;
t186 = pkin(3) * t196 + qJD(3) + t148;
t139 = t211 * t158;
t183 = t156 * t199;
t182 = qJD(4) * t203;
t181 = t158 * t188;
t178 = t158 * t185;
t177 = MDP(19) * t195;
t175 = pkin(1) * t214;
t174 = qJD(3) - t210;
t146 = pkin(2) * t179;
t172 = pkin(6) * t156 - qJ(3) * t158;
t163 = t172 * qJD(2) - t189;
t110 = t163 * qJD(1) + t146;
t145 = pkin(5) * t178;
t124 = pkin(3) * t178 + t145;
t173 = -t110 * t155 + t157 * t124;
t125 = t159 * t158 + t176;
t115 = t125 * qJD(1);
t117 = t159 * qJD(2) + t186;
t109 = t115 * t157 + t117 * t155;
t171 = t115 * t155 - t117 * t157;
t170 = t125 * t157 + t155 * t212;
t169 = -qJD(1) * t154 + t147 * t156;
t168 = -0.2e1 * qJD(2) * t123;
t164 = -qJ(3) * t190 - t189;
t114 = t164 * qJD(1) + t146;
t151 = pkin(2) * t192;
t121 = t151 + t164;
t166 = pkin(5) * t160 + qJD(1) * t121 + t114;
t165 = t122 * t156 + t159 * t190;
t134 = (-qJD(3) + t148) * qJD(2);
t135 = t148 + t174;
t162 = -t134 * t158 + (t135 * t158 + (t137 + t149) * t156) * qJD(2);
t152 = pkin(2) * t196;
t141 = t157 * t178;
t133 = qJD(2) * t139;
t132 = t149 + t150;
t131 = qJD(2) * t212;
t129 = -qJ(3) * t195 + t152;
t120 = t172 * qJD(1) + t152;
t116 = t123 * t196;
t113 = t151 + t163;
t112 = qJD(2) * t188 - t198;
t1 = [0.2e1 * t156 * MDP(4) * t178 + t213 * t214 + MDP(6) * t200 - MDP(7) * t202 + (-pkin(5) * t200 + t156 * t175) * MDP(9) + (pkin(5) * t202 + t158 * t175) * MDP(10) + t162 * MDP(11) + (t156 * t168 + t166 * t158) * MDP(12) + (-t166 * t156 + t158 * t168) * MDP(13) + (t162 * pkin(5) + t114 * t136 + t121 * t123) * MDP(14) + (-t111 * t155 * t158 + (t155 * t192 - t181) * t128) * MDP(15) + ((-t126 * t155 + t128 * t157) * t192 + (-t209 + t112 * t155 + (t126 * t157 + t128 * t155) * qJD(4)) * t158) * MDP(16) + (-t147 * t181 + t111 * t156 + (t169 * t155 + t204) * qJD(2)) * MDP(17) + (t158 * t182 - t112 * t156 + (-t126 * t158 + t169 * t157) * qJD(2)) * MDP(18) + (t147 + t196) * MDP(19) * t190 + ((-t113 * t155 + t133 * t157) * t147 - t131 * t126 + t139 * t112 + (-t122 * t191 + t173) * t156 + (-t109 * t156 - t170 * t147) * qJD(4) + (-t155 * t187 + t207 + ((-t125 * t155 + t157 * t212) * qJD(1) - t171) * qJD(2)) * t158) * MDP(20) + (t139 * t111 - t131 * t128 + (-(qJD(4) * t212 + t113) * t147 - (qJD(4) * t117 + t110) * t156) * t157 + (-(-qJD(4) * t125 + t133) * t147 + (qJD(2) * t122 + qJD(4) * t115 - t124) * t156) * t155 + (-t157 * t187 - t208 + (-t170 * qJD(1) - t109) * qJD(2)) * t158) * MDP(21); -MDP(4) * t183 + t161 * t213 + ((-t137 - t194) * t156 + (-t135 + t174) * t158) * qJD(1) * MDP(11) + (-t129 * t195 + t116) * MDP(12) + (0.2e1 * qJD(2) * qJD(3) + (t123 * t158 + t129 * t156) * qJD(1)) * MDP(13) + (-qJ(3) * t134 - qJD(3) * t137 - t123 * t129 + (-t137 * t156 + (-t135 - t210) * t158) * qJD(1) * pkin(5)) * MDP(14) + (-t128 * t203 + t209) * MDP(15) + ((-t112 - t205) * t157 + (-t111 + t206) * t155) * MDP(16) + (-t182 + t141 + (-t156 * t203 - t204) * qJD(1)) * MDP(17) + (-t147 * t188 + (-t156 * t201 + (t126 - t193) * t158) * qJD(1)) * MDP(18) - t147 * t177 + (qJ(3) * t112 + t208 - (-t120 * t155 + t132 * t157) * t147 + t186 * t126 + (t122 * t157 - t159 * t203) * qJD(4) + (t165 * t157 + t158 * t171) * qJD(1)) * MDP(20) + (qJ(3) * t111 + t207 + (t120 * t157 + t132 * t155) * t147 + t186 * t128 + (-t122 * t155 - t159 * t201) * qJD(4) + (t109 * t158 - t165 * t155) * qJD(1)) * MDP(21) + (t161 * t156 * MDP(9) + MDP(10) * t199) * pkin(1); MDP(12) * t183 + (-t153 * t161 - t160) * MDP(13) + (t116 + t145) * MDP(14) + t141 * MDP(20) + (t137 * MDP(14) - t126 * MDP(20) + (-t128 - t180) * MDP(21)) * qJD(2) + (-MDP(20) * t203 - MDP(21) * t201) * t147; t128 * t126 * MDP(15) + (-t126 ^ 2 + t128 ^ 2) * MDP(16) + (t140 + t206) * MDP(17) + (t198 + t205) * MDP(18) + qJD(2) * t177 + (t109 * t147 - t122 * t128 + t173) * MDP(20) + (-t110 * t157 + t122 * t126 - t155 * t124 - t147 * t171) * MDP(21) + (-t126 * MDP(17) - MDP(18) * t191 - t109 * MDP(20) + t171 * MDP(21)) * qJD(4);];
tauc = t1;
