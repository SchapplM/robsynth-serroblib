% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:50
% EndTime: 2019-12-31 18:02:54
% DurationCPUTime: 1.34s
% Computational Cost: add. (671->198), mult. (1391->306), div. (0->0), fcn. (796->6), ass. (0->98)
t151 = sin(pkin(8));
t159 = qJD(4) ^ 2;
t160 = qJD(1) ^ 2;
t226 = t151 * (t159 + t160);
t225 = qJ(2) * MDP(6) + MDP(5);
t155 = sin(qJ(4));
t157 = cos(qJ(4));
t201 = MDP(16) * t157;
t224 = MDP(15) * t155 + t201;
t149 = t155 ^ 2;
t223 = (-t157 ^ 2 + t149) * MDP(11);
t158 = -pkin(1) - pkin(2);
t222 = qJD(4) * pkin(7);
t139 = qJD(1) * t158 + qJD(2);
t152 = cos(pkin(8));
t200 = qJD(1) * qJ(2);
t129 = t151 * t139 + t152 * t200;
t125 = -qJD(1) * pkin(6) + t129;
t120 = qJD(3) * t155 + t125 * t157;
t187 = qJD(1) * qJD(2);
t176 = t152 * t187;
t115 = qJD(4) * t120 + t155 * t176;
t154 = sin(qJ(5));
t220 = t115 * t154;
t156 = cos(qJ(5));
t219 = t115 * t156;
t186 = qJD(1) * qJD(4);
t175 = t157 * t186;
t194 = qJD(5) * t154;
t178 = t155 * t194;
t185 = qJD(4) * qJD(5);
t206 = qJD(1) * t178 + t156 * t185;
t122 = -t156 * t175 + t206;
t218 = t122 * t154;
t217 = t122 * t157;
t143 = t154 * t185;
t193 = qJD(5) * t156;
t177 = t155 * t193;
t189 = t154 * qJD(4);
t123 = -t143 + (t157 * t189 + t177) * qJD(1);
t216 = t123 * t157;
t190 = t154 * qJD(1);
t196 = qJD(4) * t156;
t131 = t155 * t190 + t196;
t199 = qJD(1) * t157;
t140 = qJD(5) + t199;
t215 = t131 * t140;
t214 = t131 * t155;
t209 = t155 * t156;
t132 = qJD(1) * t209 - t189;
t213 = t132 * t140;
t212 = t132 * t155;
t211 = t154 * t140;
t210 = t154 * t157;
t208 = t156 * t140;
t207 = t156 * t157;
t205 = t152 * qJ(2) + t151 * t158;
t198 = qJD(2) * t152;
t134 = -pkin(6) + t205;
t197 = qJD(4) * t134;
t195 = qJD(4) * t157;
t192 = qJD(5) * t157;
t119 = qJD(3) * t157 - t125 * t155;
t116 = -qJD(4) * pkin(4) - t119;
t191 = t116 * qJD(5);
t188 = t155 * MDP(21);
t184 = 0.2e1 * t187;
t183 = 0.2e1 * t186;
t182 = t140 * t210;
t181 = t140 * t207;
t180 = t140 * t194;
t179 = t140 * t193;
t117 = t120 + t222;
t174 = t134 * t140 + t117;
t128 = t139 * t152 - t151 * t200;
t173 = -t151 * qJ(2) + t152 * t158;
t172 = t151 * t184;
t171 = 0.2e1 * t175;
t170 = t140 * t177;
t133 = pkin(3) - t173;
t169 = pkin(4) * t157 + pkin(7) * t155;
t168 = -pkin(4) * t155 + pkin(7) * t157;
t124 = qJD(1) * pkin(3) - t128;
t118 = qJD(1) * t169 + t124;
t113 = t117 * t156 + t118 * t154;
t167 = t117 * t154 - t118 * t156;
t166 = t128 * t151 - t129 * t152;
t164 = t156 * t149 * t186 + t140 * t178;
t163 = -t134 * t159 + t172;
t162 = qJD(4) * (-qJD(1) * t133 - t124 - t198);
t130 = qJD(2) * t151 + qJD(4) * t168;
t114 = qJD(4) * t119 + t157 * t176;
t161 = -qJD(4) * t116 - qJD(5) * t118 - t140 * t198 - t114;
t135 = t168 * qJD(1);
t127 = t133 + t169;
t126 = t130 * qJD(1);
t121 = t156 * t126;
t1 = [MDP(7) * t172 + 0.2e1 * MDP(8) * t176 + ((-t151 * t173 + t152 * t205) * qJD(1) - t166) * qJD(2) * MDP(9) - t183 * t223 - t159 * t157 * MDP(12) + t163 * t157 * MDP(15) + t162 * t201 + (-t122 * t209 + (t156 * t195 - t178) * t132) * MDP(17) + (-t131 * t156 - t132 * t154) * t195 * MDP(18) + (t217 + (-t181 + t212) * qJD(4) + t164) * MDP(19) + (t170 + t216 + (-t214 + (-qJD(1) * t149 + t140 * t157) * t154) * qJD(4)) * MDP(20) + (-t140 - t199) * qJD(4) * t188 + ((-t127 * t194 + t130 * t156) * t140 + (-t131 * t197 + t154 * t161 - t174 * t193 + t121) * t157) * MDP(22) + (-(t127 * t193 + t130 * t154) * t140 + (-t132 * t197 + (qJD(5) * t174 - t126) * t154 + t161 * t156) * t157) * MDP(23) + t225 * t184 + (MDP(10) * t171 + t159 * MDP(13) + t162 * MDP(15) - t163 * MDP(16) + (t218 - t123 * t156 + (t131 * t154 - t132 * t156) * qJD(5)) * MDP(18) + (-t131 * t198 - t156 * t191 - t220 - t134 * t123 + (t134 * t211 - (t127 * t156 - t134 * t210) * qJD(1) + t167) * qJD(4)) * MDP(22) + (-t132 * t198 + t154 * t191 - t219 + t134 * t122 + (t134 * t208 + (t127 * t154 + t134 * t207) * qJD(1) + t113) * qJD(4)) * MDP(23)) * t155; t166 * qJD(1) * MDP(9) + (t152 * t155 * t183 - t157 * t226) * MDP(15) + (t152 * t171 + t155 * t226) * MDP(16) + (t152 * t180 + ((t155 * t189 - t156 * t192) * t140 - t131 * t195 - t155 * t123) * t151 + (-(t151 * t156 - t152 * t210) * t140 + (-(-t151 * t210 - t152 * t156) * qJD(4) + t152 * t131) * t155) * qJD(1)) * MDP(22) + (t152 * t179 + (-(-t154 * t192 - t155 * t196) * t140 - t132 * t195 + t155 * t122) * t151 + ((t151 * t154 + t152 * t207) * t140 + ((t151 * t207 - t152 * t154) * qJD(4) + t152 * t132) * t155) * qJD(1)) * MDP(23) + (-t151 * MDP(7) - t152 * MDP(8) - t225) * t160; (-t170 + t216) * MDP(22) + (t164 - t217) * MDP(23) - t224 * t159 + ((t149 * t190 - t182 - t214) * MDP(22) + (-t181 - t212) * MDP(23)) * qJD(4); (-t132 * t208 + t218) * MDP(17) + ((t122 + t215) * t156 + (t123 + t213) * t154) * MDP(18) + t179 * MDP(19) - t180 * MDP(20) + (pkin(4) * t123 - t219 - (-t119 * t154 + t135 * t156) * t140 + t120 * t131 + (-pkin(7) * t208 + t116 * t154) * qJD(5)) * MDP(22) + (-pkin(4) * t122 + t220 + (t119 * t156 + t135 * t154) * t140 + t120 * t132 + (pkin(7) * t211 + t116 * t156) * qJD(5)) * MDP(23) + (-MDP(10) * t155 * t157 + t223) * t160 + ((t181 + (-t132 - t189) * t155) * MDP(19) + (-t182 + (t131 - t196) * t155) * MDP(20) + t140 * t188 + (-t167 * t155 + (t116 * t157 + t155 * t222) * t154) * MDP(22) + (t116 * t207 + (pkin(7) * t196 - t113) * t155) * MDP(23) + t224 * (t124 - t198)) * qJD(1); t132 * t131 * MDP(17) + (-t131 ^ 2 + t132 ^ 2) * MDP(18) + (t206 - t215) * MDP(19) + (-t143 - t213) * MDP(20) + (t113 * t140 - t114 * t154 + t116 * t132 + t121) * MDP(22) + (-t156 * t114 - t116 * t131 - t154 * t126 - t140 * t167) * MDP(23) + (-MDP(22) * t113 + MDP(23) * t167) * qJD(5) + (MDP(20) * t177 + (-t188 + (-MDP(19) * t156 + MDP(20) * t154) * t157) * qJD(4)) * qJD(1);];
tauc = t1;
