% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR9
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:58
% EndTime: 2019-12-31 20:22:09
% DurationCPUTime: 3.73s
% Computational Cost: add. (7237->392), mult. (18250->543), div. (0->0), fcn. (13277->8), ass. (0->198)
t173 = cos(qJ(2));
t238 = cos(pkin(9));
t196 = t238 * t173;
t154 = qJD(1) * t196;
t168 = sin(pkin(9));
t171 = sin(qJ(2));
t218 = qJD(1) * t171;
t130 = t168 * t218 - t154;
t172 = cos(qJ(4));
t259 = cos(qJ(5));
t206 = t259 * t172;
t169 = sin(qJ(5));
t170 = sin(qJ(4));
t229 = t169 * t170;
t183 = t206 - t229;
t265 = qJD(4) + qJD(5);
t203 = t259 * qJD(5);
t266 = qJD(4) * t259 + t203;
t252 = -t183 * t130 - t172 * t266 + t229 * t265;
t143 = t168 * t173 + t171 * t238;
t220 = qJD(1) * t143;
t107 = -qJD(2) * t172 + t170 * t220;
t109 = qJD(2) * t170 + t172 * t220;
t184 = -t107 * t169 + t109 * t259;
t48 = t107 * t259 + t109 * t169;
t257 = t48 * t184;
t207 = t259 * t170;
t146 = t169 * t172 + t207;
t100 = t265 * t146;
t240 = t146 * t130 + t100;
t159 = pkin(2) * t168 + pkin(7);
t256 = pkin(8) + t159;
t198 = qJD(4) * t256;
t227 = t170 * t130;
t211 = pkin(2) * t218;
t78 = pkin(3) * t220 + pkin(7) * t130 + t211;
t255 = -qJ(3) - pkin(6);
t152 = t255 * t173;
t148 = qJD(1) * t152;
t136 = t168 * t148;
t151 = t255 * t171;
t147 = qJD(1) * t151;
t96 = t147 * t238 + t136;
t40 = t170 * t78 + t172 * t96;
t278 = pkin(8) * t227 + t170 * t198 + t40;
t39 = -t170 * t96 + t172 * t78;
t277 = pkin(4) * t220 + t39 + (pkin(8) * t130 + t198) * t172;
t214 = qJD(1) * qJD(2);
t202 = t171 * t214;
t153 = t168 * t202;
t123 = qJD(2) * t154 - t153;
t276 = qJD(2) * qJD(4) + t123;
t217 = qJD(4) * t170;
t275 = t217 + t227;
t274 = t184 ^ 2 - t48 ^ 2;
t125 = qJD(4) + t130;
t121 = qJD(5) + t125;
t216 = qJD(4) * t172;
t208 = t170 * t276 + t220 * t216;
t215 = qJD(5) * t169;
t63 = -t172 * t276 + t220 * t217;
t21 = t107 * t203 + t109 * t215 + t169 * t208 + t259 * t63;
t273 = t121 * t48 - t21;
t164 = -pkin(2) * t173 - pkin(1);
t219 = qJD(1) * t164;
t149 = qJD(3) + t219;
t67 = t130 * pkin(3) - pkin(7) * t220 + t149;
t251 = qJD(2) * pkin(2);
t139 = t147 + t251;
t197 = t238 * t148;
t93 = t139 * t168 - t197;
t86 = qJD(2) * pkin(7) + t93;
t35 = -t170 * t86 + t172 * t67;
t30 = -pkin(8) * t109 + t35;
t25 = pkin(4) * t125 + t30;
t36 = t170 * t67 + t172 * t86;
t31 = -pkin(8) * t107 + t36;
t132 = t143 * qJD(2);
t122 = qJD(1) * t132;
t199 = qJD(2) * t255;
t128 = t173 * qJD(3) + t171 * t199;
t115 = t128 * qJD(1);
t178 = -t171 * qJD(3) + t173 * t199;
t116 = t178 * qJD(1);
t62 = t115 * t238 + t116 * t168;
t157 = pkin(2) * t202;
t66 = pkin(3) * t122 - pkin(7) * t123 + t157;
t17 = -qJD(4) * t36 - t170 * t62 + t172 * t66;
t6 = t122 * pkin(4) + t63 * pkin(8) + t17;
t16 = t170 * t66 + t172 * t62 + t216 * t67 - t217 * t86;
t9 = -pkin(8) * t208 + t16;
t177 = -t169 * t6 - t203 * t25 + t215 * t31 - t259 * t9;
t92 = t139 * t238 + t136;
t85 = -qJD(2) * pkin(3) - t92;
t46 = pkin(4) * t107 + t85;
t272 = t46 * t48 + t177;
t270 = -0.2e1 * t214;
t269 = -t125 * t35 + t16;
t268 = t36 * t125 + t17;
t194 = t125 * t170;
t267 = t109 * t194;
t210 = t259 * t31;
t8 = t169 * t25 + t210;
t2 = -qJD(5) * t8 - t169 * t9 + t259 * t6;
t264 = -t184 * t46 + t2;
t22 = qJD(5) * t184 - t169 * t63 + t208 * t259;
t263 = t121 * t184 - t22;
t262 = -t121 * t252 + t146 * t122;
t261 = -t183 * t21 - t184 * t240;
t260 = t220 ^ 2;
t258 = pkin(2) * t171;
t140 = t256 * t170;
t141 = t256 * t172;
t89 = -t140 * t169 + t141 * t259;
t254 = t89 * qJD(5) - t169 * t278 + t259 * t277;
t88 = -t140 * t259 - t141 * t169;
t253 = -t88 * qJD(5) + t169 * t277 + t259 * t278;
t180 = -t168 * t171 + t196;
t91 = -pkin(3) * t180 - pkin(7) * t143 + t164;
t106 = t151 * t168 - t152 * t238;
t98 = t172 * t106;
t45 = t170 * t91 + t98;
t250 = t220 * t48;
t248 = t169 * t31;
t245 = t184 * t220;
t105 = -t151 * t238 - t152 * t168;
t61 = t115 * t168 - t116 * t238;
t244 = t61 * t105;
t243 = t61 * t170;
t242 = t61 * t172;
t241 = t63 * t170;
t54 = t170 * t208;
t239 = -t107 * t216 - t54;
t237 = t107 * t130;
t236 = t109 * t107;
t235 = t109 * t220;
t94 = t122 * t180;
t234 = t220 * t107;
t233 = t220 * t130;
t232 = t143 * t170;
t231 = t143 * t172;
t228 = t170 * t122;
t135 = t180 * qJD(2);
t226 = t170 * t135;
t112 = t172 * t122;
t225 = t172 * t135;
t175 = qJD(1) ^ 2;
t224 = t173 * t175;
t174 = qJD(2) ^ 2;
t223 = t174 * t171;
t222 = t174 * t173;
t221 = t171 ^ 2 - t173 ^ 2;
t212 = t171 * t251;
t209 = t171 * t224;
t205 = t143 * t217;
t77 = t128 * t238 + t168 * t178;
t79 = pkin(3) * t132 - pkin(7) * t135 + t212;
t200 = -t170 * t77 + t172 * t79;
t44 = -t106 * t170 + t172 * t91;
t195 = pkin(1) * t270;
t76 = t168 * t128 - t178 * t238;
t95 = t147 * t168 - t197;
t193 = t125 * t172;
t192 = 0.2e1 * t220;
t191 = -t146 * t22 + t252 * t48;
t190 = -t121 * t240 + t183 * t122;
t189 = t173 * t202;
t188 = pkin(4) * t275 - t95;
t160 = -pkin(2) * t238 - pkin(3);
t187 = t208 * t172;
t186 = -t170 * t36 - t172 * t35;
t185 = -t125 * t275 + t112;
t32 = -pkin(4) * t180 - pkin(8) * t231 + t44;
t38 = -pkin(8) * t232 + t45;
t18 = -t169 * t38 + t259 * t32;
t19 = t169 * t32 + t259 * t38;
t182 = t143 * t216 + t226;
t181 = -t205 + t225;
t23 = -t106 * t217 + t170 * t79 + t172 * t77 + t216 * t91;
t179 = -t122 * t159 + t125 * t85;
t150 = -pkin(4) * t172 + t160;
t129 = t130 ^ 2;
t82 = t183 * t143;
t81 = t146 * t143;
t75 = pkin(4) * t232 + t105;
t41 = pkin(4) * t182 + t76;
t33 = pkin(4) * t208 + t61;
t27 = t135 * t207 - t169 * t205 - t215 * t232 + (t135 * t169 + t143 * t266) * t172;
t26 = t100 * t143 - t135 * t206 + t169 * t226;
t24 = -qJD(4) * t45 + t200;
t15 = -pkin(8) * t182 + t23;
t14 = -pkin(8) * t225 + t132 * pkin(4) + (-t98 + (pkin(8) * t143 - t91) * t170) * qJD(4) + t200;
t11 = t259 * t30 - t248;
t10 = -t169 * t30 - t210;
t7 = t25 * t259 - t248;
t4 = -qJD(5) * t19 + t14 * t259 - t169 * t15;
t3 = qJD(5) * t18 + t169 * t14 + t15 * t259;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t189, t221 * t270, t222, -0.2e1 * t189, -t223, 0, -pkin(6) * t222 + t171 * t195, pkin(6) * t223 + t173 * t195, 0, 0, t123 * t143 + t135 * t220, -t122 * t143 + t123 * t180 - t130 * t135 - t132 * t220, t135 * qJD(2), t130 * t132 - t94, -t132 * qJD(2), 0, t164 * t122 + t149 * t132 + (-t76 + (-qJD(1) * t180 + t130) * t258) * qJD(2), t164 * t123 + t149 * t135 + (t192 * t258 - t77) * qJD(2), t105 * t123 - t106 * t122 - t130 * t77 - t132 * t93 - t135 * t92 + t143 * t61 + t180 * t62 + t220 * t76, t244 + t62 * t106 - t92 * t76 + t93 * t77 + (t149 + t219) * t212, t109 * t181 - t231 * t63, (-t107 * t172 - t109 * t170) * t135 + (-t187 + t241 + (t107 * t170 - t109 * t172) * qJD(4)) * t143, t109 * t132 + t112 * t143 + t125 * t181 + t180 * t63, t107 * t182 + t143 * t54, -t107 * t132 - t125 * t182 - t143 * t228 + t180 * t208, t125 * t132 - t94, t24 * t125 + t44 * t122 - t17 * t180 + t35 * t132 + t76 * t107 + t105 * t208 + t85 * t226 + (t216 * t85 + t243) * t143, t85 * t225 - t105 * t63 + t76 * t109 - t45 * t122 - t23 * t125 - t36 * t132 + t16 * t180 + (-t217 * t85 + t242) * t143, -t23 * t107 - t45 * t208 - t24 * t109 + t44 * t63 + t186 * t135 + (-t16 * t170 - t17 * t172 + (t170 * t35 - t172 * t36) * qJD(4)) * t143, t16 * t45 + t17 * t44 + t23 * t36 + t24 * t35 + t76 * t85 + t244, -t184 * t26 - t21 * t82, -t184 * t27 + t21 * t81 - t22 * t82 + t26 * t48, -t121 * t26 + t122 * t82 + t132 * t184 + t180 * t21, t22 * t81 + t27 * t48, -t121 * t27 - t122 * t81 - t132 * t48 + t180 * t22, t121 * t132 - t94, t121 * t4 + t122 * t18 + t132 * t7 - t180 * t2 + t22 * t75 + t27 * t46 + t33 * t81 + t41 * t48, -t121 * t3 - t122 * t19 - t132 * t8 - t177 * t180 + t184 * t41 - t21 * t75 - t26 * t46 + t33 * t82, t177 * t81 + t18 * t21 - t184 * t4 - t19 * t22 - t2 * t82 + t26 * t7 - t27 * t8 - t3 * t48, -t177 * t19 + t18 * t2 + t3 * t8 + t33 * t75 + t4 * t7 + t41 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t221 * t175, 0, t209, 0, 0, t175 * pkin(1) * t171, pkin(1) * t224, 0, 0, t233, -t129 + t260, -t153 + (t154 + t130) * qJD(2), -t233, 0, 0, qJD(2) * t95 - t130 * t211 - t149 * t220 - t61, qJD(2) * t96 + t130 * t149 - t211 * t220 - t62, (t93 - t95) * t220 + (-t92 + t96) * t130 + (-t122 * t168 - t123 * t238) * pkin(2), t92 * t95 - t93 * t96 + (-t149 * t218 + t168 * t62 - t238 * t61) * pkin(2), t109 * t193 - t241, (-t63 - t237) * t172 - t267 + t239, t125 * t193 + t228 - t235, t107 * t194 - t187, t185 + t234, -t125 * t220, t160 * t208 - t242 - t35 * t220 - t95 * t107 + (-t159 * t216 - t39) * t125 + t179 * t170, -t95 * t109 + t36 * t220 - t160 * t63 + t243 + (t159 * t217 + t40) * t125 + t179 * t172, t40 * t107 + t39 * t109 + (-t159 * t208 + t16 - t35 * t130 + (t109 * t159 - t35) * qJD(4)) * t172 + (-t36 * t130 - t159 * t63 - t17 + (t107 * t159 - t36) * qJD(4)) * t170, t61 * t160 - t35 * t39 - t36 * t40 - t85 * t95 + (qJD(4) * t186 + t16 * t172 - t17 * t170) * t159, -t21 * t146 - t184 * t252, t191 + t261, -t245 + t262, -t183 * t22 + t240 * t48, t190 + t250, -t121 * t220, -t121 * t254 + t88 * t122 + t150 * t22 - t183 * t33 + t188 * t48 - t220 * t7 + t240 * t46, t121 * t253 - t89 * t122 + t33 * t146 - t150 * t21 + t184 * t188 + t220 * t8 - t252 * t46, -t2 * t146 - t177 * t183 + t184 * t254 + t88 * t21 - t89 * t22 - t240 * t8 + t252 * t7 + t253 * t48, t33 * t150 - t177 * t89 + t188 * t46 + t2 * t88 - t253 * t8 - t254 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192 * qJD(2), -t153 + (t154 - t130) * qJD(2), -t129 - t260, t130 * t93 + t220 * t92 + t157, 0, 0, 0, 0, 0, 0, t185 - t234, -t125 ^ 2 * t172 - t228 - t235, (t63 - t237) * t172 + t267 + t239, t170 * t269 + t172 * t268 - t220 * t85, 0, 0, 0, 0, 0, 0, t190 - t250, -t245 - t262, t191 - t261, -t146 * t177 + t183 * t2 - t220 * t46 - t240 * t7 - t252 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, -t107 ^ 2 + t109 ^ 2, t107 * t125 - t63, -t236, t109 * t125 - t208, t122, -t85 * t109 + t268, t107 * t85 - t269, 0, 0, t257, t274, t273, -t257, t263, t122, -t10 * t121 + (-t109 * t48 - t121 * t215 + t122 * t259) * pkin(4) + t264, t11 * t121 + (-t109 * t184 - t121 * t203 - t122 * t169) * pkin(4) + t272, t10 * t184 + t11 * t48 + t8 * t184 - t7 * t48 + (t259 * t21 - t169 * t22 + (t169 * t184 - t259 * t48) * qJD(5)) * pkin(4), -t7 * t10 - t8 * t11 + (t259 * t2 - t177 * t169 - t109 * t46 + (-t169 * t7 + t259 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, t274, t273, -t257, t263, t122, t8 * t121 + t264, t121 * t7 + t272, 0, 0;];
tauc_reg = t1;
