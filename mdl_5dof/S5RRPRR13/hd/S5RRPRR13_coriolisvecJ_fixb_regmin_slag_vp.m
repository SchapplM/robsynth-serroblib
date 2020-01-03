% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR13_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:52
% EndTime: 2019-12-31 20:34:01
% DurationCPUTime: 3.18s
% Computational Cost: add. (3259->330), mult. (8418->484), div. (0->0), fcn. (6337->8), ass. (0->180)
t175 = cos(qJ(2));
t218 = t175 * qJD(1);
t157 = -qJD(4) + t218;
t152 = -qJD(5) + t157;
t170 = sin(qJ(5));
t173 = cos(qJ(5));
t168 = sin(pkin(9));
t172 = sin(qJ(2));
t229 = qJD(1) * t172;
t212 = t168 * t229;
t169 = cos(pkin(9));
t219 = t169 * qJD(2);
t127 = t212 - t219;
t211 = t169 * t229;
t220 = t168 * qJD(2);
t129 = t211 + t220;
t171 = sin(qJ(4));
t174 = cos(qJ(4));
t75 = t171 * t127 - t174 * t129;
t76 = t174 * t127 + t171 * t129;
t253 = t170 * t75 - t173 * t76;
t264 = t253 * t152;
t194 = pkin(2) * t172 - qJ(3) * t175;
t137 = t194 * qJD(1);
t118 = t168 * t137;
t240 = t169 * t172;
t241 = t168 * t175;
t186 = -pkin(6) * t240 - pkin(7) * t241;
t82 = qJD(1) * t186 + t118;
t263 = qJD(3) * t169 - t82;
t192 = t170 * t76 + t173 * t75;
t262 = t192 * t253;
t261 = t75 * t157;
t260 = t76 * t157;
t259 = t152 * t192;
t237 = t174 * t169;
t238 = t171 * t168;
t134 = -t237 + t238;
t184 = t134 * t175;
t233 = qJD(1) * t184 - t134 * qJD(4);
t135 = t174 * t168 + t171 * t169;
t185 = t135 * t175;
t232 = -qJD(1) * t185 + t135 * qJD(4);
t258 = t192 ^ 2 - t253 ^ 2;
t214 = pkin(3) * t218;
t142 = -t175 * pkin(2) - t172 * qJ(3) - pkin(1);
t120 = t142 * qJD(1);
t162 = pkin(6) * t218;
t147 = qJD(2) * qJ(3) + t162;
t84 = t169 * t120 - t168 * t147;
t45 = -t129 * pkin(7) - t214 + t84;
t85 = t168 * t120 + t169 * t147;
t49 = -t127 * pkin(7) + t85;
t18 = t171 * t45 + t174 * t49;
t13 = -t76 * pkin(8) + t18;
t222 = qJD(5) * t170;
t11 = t13 * t222;
t161 = pkin(6) * t229;
t244 = qJD(2) * pkin(2);
t202 = qJD(3) - t244;
t141 = t161 + t202;
t91 = t127 * pkin(3) + t141;
t35 = t76 * pkin(4) + t91;
t257 = -t253 * t35 + t11;
t221 = qJD(5) * t173;
t217 = qJD(1) * qJD(2);
t209 = t175 * t217;
t198 = t168 * t209;
t223 = qJD(4) * t174;
t38 = -t127 * t223 + t209 * t237 + (-qJD(4) * t129 - t198) * t171;
t182 = qJD(2) * t185;
t251 = qJD(4) * t75;
t39 = qJD(1) * t182 - t251;
t7 = -t170 * t39 + t173 * t38 - t76 * t221 + t222 * t75;
t256 = t7 + t264;
t17 = -t171 * t49 + t174 * t45;
t12 = pkin(8) * t75 + t17;
t10 = -t157 * pkin(4) + t12;
t243 = t173 * t13;
t2 = t170 * t10 + t243;
t160 = t172 * t217;
t239 = t169 * t175;
t189 = pkin(3) * t172 - pkin(7) * t239;
t183 = t189 * qJD(2);
t113 = qJD(2) * t194 - t172 * qJD(3);
t102 = t113 * qJD(1);
t140 = (qJD(3) - t161) * qJD(2);
t63 = t169 * t102 - t168 * t140;
t42 = qJD(1) * t183 + t63;
t64 = t168 * t102 + t169 * t140;
t46 = -pkin(7) * t198 + t64;
t206 = -t171 * t46 + t174 * t42;
t179 = -qJD(4) * t18 + t206;
t4 = pkin(4) * t160 - t38 * pkin(8) + t179;
t225 = qJD(4) * t171;
t188 = t171 * t42 + t174 * t46 + t45 * t223 - t49 * t225;
t5 = -t39 * pkin(8) + t188;
t213 = -t170 * t5 + t173 * t4;
t255 = -qJD(5) * t2 + t35 * t192 + t213;
t178 = qJD(5) * t192 - t170 * t38 - t173 * t39;
t254 = t178 + t259;
t252 = -0.2e1 * t217;
t249 = pkin(7) + qJ(3);
t144 = t249 * t168;
t145 = t249 * t169;
t231 = -t171 * t144 + t174 * t145;
t190 = qJD(3) * t168 + qJD(4) * t145;
t92 = pkin(6) * t212 + t169 * t137;
t65 = qJD(1) * t189 + t92;
t250 = -t144 * t223 + t263 * t174 + (-t190 - t65) * t171;
t80 = t173 * t134 + t170 * t135;
t248 = -qJD(5) * t80 - t232 * t170 + t233 * t173;
t81 = -t170 * t134 + t173 * t135;
t247 = qJD(5) * t81 + t233 * t170 + t232 * t173;
t126 = t169 * t142;
t83 = -pkin(7) * t240 + t126 + (-pkin(6) * t168 - pkin(3)) * t175;
t242 = t168 * t172;
t155 = pkin(6) * t239;
t98 = t168 * t142 + t155;
t90 = -pkin(7) * t242 + t98;
t245 = t171 * t83 + t174 * t90;
t177 = qJD(1) ^ 2;
t236 = t175 * t177;
t176 = qJD(2) ^ 2;
t235 = t176 * t172;
t234 = t176 * t175;
t228 = qJD(2) * t172;
t215 = pkin(6) * t228;
t88 = t169 * t113 + t168 * t215;
t156 = pkin(6) * t209;
t112 = pkin(3) * t198 + t156;
t121 = t168 * t214 + t162;
t227 = qJD(2) * t175;
t163 = pkin(6) * t227;
t122 = t175 * pkin(3) * t220 + t163;
t138 = pkin(3) * t242 + t172 * pkin(6);
t230 = t172 ^ 2 - t175 ^ 2;
t224 = qJD(4) * t172;
t216 = pkin(6) * t241;
t159 = -t169 * pkin(3) - pkin(2);
t210 = t232 * pkin(4) - t121;
t208 = qJD(5) * t10 + t5;
t55 = t183 + t88;
t104 = t168 * t113;
t67 = qJD(2) * t186 + t104;
t205 = -t171 * t67 + t174 * t55;
t204 = -t171 * t90 + t174 * t83;
t203 = pkin(1) * t252;
t201 = -t174 * t144 - t171 * t145;
t200 = t127 + t219;
t199 = -t129 + t220;
t57 = t174 * t65;
t59 = -t134 * pkin(8) + t231;
t197 = pkin(4) * t229 + t233 * pkin(8) + t135 * qJD(3) + t231 * qJD(4) + qJD(5) * t59 - t171 * t82 + t57;
t58 = -t135 * pkin(8) + t201;
t196 = -t232 * pkin(8) + qJD(5) * t58 + t250;
t195 = -t141 + t202;
t110 = t134 * t172;
t20 = -t175 * pkin(4) + t110 * pkin(8) + t204;
t109 = t135 * t172;
t21 = -t109 * pkin(8) + t245;
t193 = t170 * t20 + t173 * t21;
t53 = t173 * t109 - t170 * t110;
t54 = -t170 * t109 - t173 * t110;
t187 = t171 * t55 + t174 * t67 + t83 * t223 - t90 * t225;
t103 = t134 * pkin(4) + t159;
t97 = t126 - t216;
t93 = -pkin(6) * t211 + t118;
t89 = -t169 * t215 + t104;
t87 = t109 * pkin(4) + t138;
t61 = t223 * t240 - t224 * t238 + t182;
t60 = -qJD(2) * t184 - t135 * t224;
t40 = t61 * pkin(4) + t122;
t22 = t39 * pkin(4) + t112;
t15 = qJD(5) * t54 + t170 * t60 + t173 * t61;
t14 = -qJD(5) * t53 - t170 * t61 + t173 * t60;
t9 = -t61 * pkin(8) + t187;
t6 = pkin(4) * t228 - t60 * pkin(8) - qJD(4) * t245 + t205;
t1 = t173 * t10 - t170 * t13;
t3 = [0, 0, 0, 0.2e1 * t175 * t160, t230 * t252, t234, -t235, 0, -pkin(6) * t234 + t172 * t203, pkin(6) * t235 + t175 * t203, (-qJD(1) * t88 - t63) * t175 + ((pkin(6) * t127 + t141 * t168) * t175 + (t84 + (t97 + 0.2e1 * t216) * qJD(1)) * t172) * qJD(2), (qJD(1) * t89 + t64) * t175 + ((pkin(6) * t129 + t141 * t169) * t175 + (-t85 + (-t98 + 0.2e1 * t155) * qJD(1)) * t172) * qJD(2), -t89 * t127 - t88 * t129 + (-t168 * t64 - t169 * t63) * t172 + (-t168 * t85 - t169 * t84 + (-t168 * t98 - t169 * t97) * qJD(1)) * t227, t63 * t97 + t64 * t98 + t84 * t88 + t85 * t89 + (t141 + t161) * t163, -t38 * t110 - t60 * t75, -t38 * t109 + t110 * t39 - t60 * t76 + t61 * t75, -t60 * t157 - t38 * t175 + (-qJD(1) * t110 - t75) * t228, t61 * t157 + t39 * t175 + (-qJD(1) * t109 - t76) * t228, (-t157 - t218) * t228, -t205 * t157 - t206 * t175 + t122 * t76 + t138 * t39 + t112 * t109 + t91 * t61 + (t157 * t245 + t175 * t18) * qJD(4) + (qJD(1) * t204 + t17) * t228, t187 * t157 + t188 * t175 - t122 * t75 + t138 * t38 - t112 * t110 + t91 * t60 + (-t245 * qJD(1) - t18) * t228, -t14 * t192 + t7 * t54, t14 * t253 + t15 * t192 + t178 * t54 - t7 * t53, -t14 * t152 - t7 * t175 + (qJD(1) * t54 - t192) * t228, t15 * t152 - t178 * t175 + (-qJD(1) * t53 + t253) * t228, (-t152 - t218) * t228, -(-t170 * t9 + t173 * t6) * t152 - t213 * t175 - t40 * t253 - t87 * t178 + t22 * t53 + t35 * t15 + (t152 * t193 + t175 * t2) * qJD(5) + ((-t170 * t21 + t173 * t20) * qJD(1) + t1) * t228, -t11 * t175 + t35 * t14 + t22 * t54 - t40 * t192 + t87 * t7 + ((-qJD(5) * t21 + t6) * t152 + t4 * t175) * t170 + ((qJD(5) * t20 + t9) * t152 + t208 * t175) * t173 + (-qJD(1) * t193 - t2) * t228; 0, 0, 0, -t172 * t236, t230 * t177, 0, 0, 0, t177 * pkin(1) * t172, pkin(1) * t236, ((-qJ(3) * t220 - t84) * t172 + (-pkin(6) * t200 + t168 * t195 + t92) * t175) * qJD(1), ((-qJ(3) * t219 + t85) * t172 + (pkin(6) * t199 + t169 * t195 - t93) * t175) * qJD(1), t93 * t127 + t92 * t129 + (-qJD(3) * t127 + t84 * t218 + t64) * t169 + (qJD(3) * t129 + t85 * t218 - t63) * t168, -t84 * t92 - t85 * t93 + (-t168 * t84 + t169 * t85) * qJD(3) + (-t63 * t168 + t64 * t169) * qJ(3) + (-t141 - t244) * t162, t38 * t135 - t233 * t75, -t38 * t134 - t135 * t39 + t232 * t75 - t233 * t76, -t233 * t157 + (qJD(2) * t135 + t75) * t229, t232 * t157 + (-qJD(2) * t134 + t76) * t229, t157 * t229, t112 * t134 - t121 * t76 + t159 * t39 + t232 * t91 + (t57 + t190 * t174 + (-qJD(4) * t144 + t263) * t171) * t157 + (qJD(2) * t201 - t17) * t229, t112 * t135 + t121 * t75 + t159 * t38 + t233 * t91 + t250 * t157 + (-t231 * qJD(2) + t18) * t229, -t192 * t248 + t7 * t81, t178 * t81 + t192 * t247 + t248 * t253 - t7 * t80, -t248 * t152 + (qJD(2) * t81 + t192) * t229, t247 * t152 + (-qJD(2) * t80 - t253) * t229, t152 * t229, -t103 * t178 + t22 * t80 + t247 * t35 - t210 * t253 + (t170 * t196 + t173 * t197) * t152 + ((-t170 * t59 + t173 * t58) * qJD(2) - t1) * t229, t103 * t7 + t22 * t81 + t248 * t35 - t210 * t192 + (-t170 * t197 + t173 * t196) * t152 + (-(t170 * t58 + t173 * t59) * qJD(2) + t2) * t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199 * t218, t200 * t218, -t127 ^ 2 - t129 ^ 2, t85 * t127 + t84 * t129 + t156, 0, 0, 0, 0, 0, t39 + t261, t38 + t260, 0, 0, 0, 0, 0, -t178 + t259, t7 - t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75 * t76, t75 ^ 2 - t76 ^ 2, t38 - t260, -t135 * t209 + t251 + t261, t160, -t18 * t157 + t75 * t91 + t179, -t17 * t157 + t91 * t76 - t188, t262, t258, t256, t254, t160, (-t170 * t12 - t243) * t152 + (t152 * t222 + t160 * t173 - t253 * t75) * pkin(4) + t255, (t13 * t152 - t4) * t170 + (-t12 * t152 - t208) * t173 + (t152 * t221 - t160 * t170 - t192 * t75) * pkin(4) + t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, t258, t256, t254, t160, -t2 * t152 + t255, -t1 * t152 - t170 * t4 - t173 * t208 + t257;];
tauc_reg = t3;
