% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:07:43
% EndTime: 2019-03-08 22:07:58
% DurationCPUTime: 5.60s
% Computational Cost: add. (5783->424), mult. (17613->629), div. (0->0), fcn. (14949->14), ass. (0->227)
t166 = cos(pkin(7));
t174 = cos(qJ(3));
t175 = cos(qJ(2));
t262 = t174 * t175;
t170 = sin(qJ(3));
t171 = sin(qJ(2));
t266 = t170 * t171;
t188 = -t166 * t266 + t262;
t252 = qJD(3) * t174;
t228 = t166 * t252;
t164 = sin(pkin(6));
t259 = qJD(1) * t164;
t308 = pkin(2) * t228 - t188 * t259;
t163 = sin(pkin(7));
t294 = pkin(9) + qJ(4);
t226 = t294 * t170;
t250 = qJD(4) * t174;
t307 = (-qJD(3) * t226 + t250) * t163 + t308;
t264 = t171 * t174;
t265 = t170 * t175;
t190 = -t166 * t264 - t265;
t123 = t190 * t259;
t270 = t166 * t170;
t273 = t163 * t174;
t128 = pkin(2) * t270 + t294 * t273;
t251 = qJD(4) * t170;
t306 = -t128 * qJD(3) - t163 * t251 - t123;
t162 = sin(pkin(13));
t165 = cos(pkin(13));
t290 = t306 * t162 + t307 * t165;
t195 = t162 * t174 + t165 * t170;
t136 = t195 * t163;
t131 = qJD(3) * t136;
t254 = qJD(3) * t163;
t271 = t165 * t174;
t132 = (-t162 * t170 + t271) * t254;
t253 = qJD(3) * t170;
t230 = t163 * t253;
t233 = t171 * t259;
t305 = -pkin(3) * t230 - pkin(4) * t131 + pkin(10) * t132 + t163 * t233;
t257 = qJD(2) * t163;
t133 = t195 * t257;
t173 = cos(qJ(5));
t243 = t166 * qJD(2);
t216 = qJD(3) + t243;
t151 = t173 * t216;
t169 = sin(qJ(5));
t98 = t133 * t169 - t151;
t97 = qJD(6) + t98;
t231 = t170 * t257;
t238 = t163 * t271;
t130 = qJD(2) * t238 - t162 * t231;
t126 = qJD(5) - t130;
t291 = t307 * t162 - t306 * t165;
t168 = sin(qJ(6));
t172 = cos(qJ(6));
t182 = -t173 * t133 - t169 * t216;
t63 = -t172 * t126 - t168 * t182;
t304 = t126 * t63;
t246 = qJD(5) * t173;
t248 = qJD(5) * t169;
t115 = (pkin(2) * t174 + pkin(3)) * t166 - t163 * t226;
t82 = t162 * t115 + t165 * t128;
t69 = pkin(10) * t166 + t82;
t274 = t163 * t170;
t135 = t162 * t274 - t238;
t235 = -pkin(3) * t174 - pkin(2);
t92 = pkin(4) * t135 - pkin(10) * t136 + t235 * t163;
t303 = t305 * t169 - t173 * t290 - t92 * t246 + t69 * t248;
t302 = -t169 * t290 - t305 * t173;
t301 = t170 * t174;
t232 = t175 * t259;
t289 = qJD(2) * pkin(2);
t148 = t232 + t289;
t167 = cos(pkin(6));
t258 = qJD(1) * t167;
t234 = t163 * t258;
t185 = -t148 * t166 - t234;
t141 = pkin(9) * t257 + t233;
t263 = t174 * t141;
t300 = (t185 * t170 - t263) * qJD(3) + ((-qJ(4) * t252 - t251) * t163 + t123) * qJD(2);
t121 = qJD(2) * t131;
t255 = qJD(2) * t174;
t91 = t148 * t270 + t263 + (qJ(4) * t255 + t170 * t258) * t163;
t285 = t165 * t91;
t149 = t174 * t234;
t269 = t166 * t174;
t261 = t148 * t269 + t149;
t90 = (-qJ(4) * t257 - t141) * t170 + t261;
t67 = t216 * pkin(3) + t90;
t34 = t162 * t67 + t285;
t32 = t216 * pkin(10) + t34;
t155 = t166 * t258;
t104 = qJD(4) + t155 + (-pkin(3) * t255 - t148) * t163;
t46 = -pkin(4) * t130 - pkin(10) * t133 + t104;
t16 = t169 * t46 + t173 * t32;
t256 = qJD(2) * t164;
t225 = qJD(1) * t256;
t209 = t175 * t225;
t236 = qJD(3) * t149 + t148 * t228 + t174 * t209;
t42 = -t141 * t253 + (-t233 * t270 + (-qJ(4) * t253 + t250) * t163) * qJD(2) + t236;
t22 = t300 * t162 + t165 * t42;
t242 = qJD(2) * qJD(3);
t224 = t163 * t242;
t207 = t174 * t224;
t208 = t170 * t224;
t122 = -t162 * t208 + t165 * t207;
t210 = t171 * t225;
t127 = pkin(3) * t208 + t163 * t210;
t53 = pkin(4) * t121 - pkin(10) * t122 + t127;
t222 = t169 * t22 - t173 * t53;
t4 = -pkin(5) * t121 + t16 * qJD(5) + t222;
t299 = (-pkin(5) * t182 + t97 * pkin(11)) * t97 + t4;
t61 = qJD(5) * t151 + t173 * t122 - t133 * t248;
t65 = t126 * t168 - t172 * t182;
t26 = t65 * qJD(6) - t172 * t121 + t168 * t61;
t62 = -t182 * qJD(5) + t169 * t122;
t12 = pkin(11) * t126 + t16;
t73 = t162 * t91;
t33 = t165 * t67 - t73;
t31 = -t216 * pkin(4) - t33;
t17 = t98 * pkin(5) + pkin(11) * t182 + t31;
t202 = t12 * t168 - t17 * t172;
t187 = t169 * t53 + t173 * t22 + t46 * t246 - t32 * t248;
t3 = pkin(11) * t121 + t187;
t21 = t162 * t42 - t165 * t300;
t8 = pkin(5) * t62 - pkin(11) * t61 + t21;
t1 = -t202 * qJD(6) + t168 * t8 + t172 * t3;
t176 = qJD(2) ^ 2;
t298 = t63 * t97;
t297 = t65 * t97;
t276 = t130 * t173;
t93 = -t172 * t133 + t168 * t276;
t296 = t93 * t97;
t94 = t133 * t168 + t172 * t276;
t295 = t94 * t97;
t198 = t169 * t92 + t173 * t69;
t293 = -pkin(5) * t131 + t198 * qJD(5) - t302;
t40 = t165 * t90 - t73;
t86 = pkin(3) * t231 + pkin(4) * t133 - pkin(10) * t130;
t292 = t169 * t86 + t173 * t40;
t288 = t126 * t98;
t287 = t130 * t65;
t286 = t133 * t98;
t244 = qJD(6) * t172;
t245 = qJD(6) * t168;
t25 = t168 * t121 + t126 * t244 + t172 * t61 + t182 * t245;
t284 = t168 * t25;
t283 = t168 * t62;
t282 = t172 * t62;
t220 = t172 * t97;
t39 = t162 * t90 + t285;
t281 = -t39 + t126 * (pkin(5) * t169 - pkin(11) * t173);
t280 = t182 * t126;
t279 = t182 * t133;
t125 = -t148 * t163 + t155;
t278 = t125 * t163;
t277 = t126 * t169;
t159 = t163 ^ 2;
t275 = t159 * t176;
t272 = t164 * t176;
t268 = t169 * t121;
t260 = t170 ^ 2 - t174 ^ 2;
t249 = qJD(5) * t168;
t247 = qJD(5) * t172;
t240 = t97 * t249;
t239 = t97 * t247;
t237 = t171 * t272;
t158 = -pkin(3) * t165 - pkin(4);
t229 = t163 * t252;
t227 = t163 * t166 * t176;
t140 = -pkin(5) * t173 - pkin(11) * t169 + t158;
t219 = pkin(11) * t133 - qJD(6) * t140 + t292;
t81 = t115 * t165 - t162 * t128;
t218 = t126 * t173;
t217 = 0.2e1 * t159 * t242;
t215 = qJD(3) + 0.2e1 * t243;
t214 = t159 * t237;
t212 = t163 * t171 * t256;
t107 = t136 * t169 - t173 * t166;
t108 = t136 * t173 + t166 * t169;
t68 = -pkin(4) * t166 - t81;
t37 = pkin(5) * t107 - pkin(11) * t108 + t68;
t211 = -pkin(11) * t131 - qJD(6) * t37 + t303;
t30 = pkin(11) * t135 + t198;
t79 = -t107 * qJD(5) + t132 * t173;
t80 = t108 * qJD(5) + t132 * t169;
t206 = -pkin(5) * t80 + pkin(11) * t79 + qJD(6) * t30 - t291;
t203 = -t141 + t233;
t6 = t12 * t172 + t168 * t17;
t139 = -t163 * t164 * t175 + t166 * t167;
t189 = t166 * t265 + t264;
t105 = t189 * t164 + t167 * t274;
t191 = t166 * t262 - t266;
t179 = t191 * t164 + t167 * t273;
t59 = t165 * t105 + t162 * t179;
t44 = t139 * t169 + t173 * t59;
t58 = t105 * t162 - t165 * t179;
t201 = t168 * t58 + t172 * t44;
t200 = -t168 * t44 + t172 * t58;
t15 = -t169 * t32 + t173 * t46;
t197 = -t169 * t69 + t173 * t92;
t196 = t139 * t173 - t169 * t59;
t89 = t108 * t172 + t135 * t168;
t88 = t108 * t168 - t172 * t135;
t194 = t173 * t121 - t126 * t248 + t130 * t277;
t193 = -t97 * t244 - t283;
t192 = t97 * t245 - t282;
t157 = pkin(3) * t162 + pkin(10);
t184 = -t157 * t121 + t126 * t31;
t11 = -pkin(5) * t126 - t15;
t183 = -pkin(11) * t62 + (t11 + t15) * t97;
t181 = -qJD(3) * t141 - t166 * t210;
t2 = -t6 * qJD(6) - t168 * t3 + t172 * t8;
t180 = t185 * t166 + t278;
t78 = t167 * t229 + (t188 * qJD(2) + t191 * qJD(3)) * t164;
t77 = -t167 * t230 + (t190 * qJD(2) - t189 * qJD(3)) * t164;
t56 = t65 * t248;
t36 = t162 * t77 + t165 * t78;
t35 = t162 * t78 - t165 * t77;
t29 = -pkin(5) * t135 - t197;
t28 = t89 * qJD(6) - t172 * t131 + t168 * t79;
t27 = -t88 * qJD(6) + t131 * t168 + t172 * t79;
t18 = -pkin(5) * t133 + t169 * t40 - t173 * t86;
t14 = t44 * qJD(5) + t169 * t36 - t173 * t212;
t13 = t196 * qJD(5) + t169 * t212 + t173 * t36;
t5 = [0, 0, -t237, -t175 * t272, 0, 0, 0, 0, 0, t139 * t208 - t174 * t214 + t216 * t77, t139 * t207 + t170 * t214 - t216 * t78, -t121 * t59 + t122 * t58 + t130 * t36 + t133 * t35, t104 * t212 + t127 * t139 + t21 * t58 + t22 * t59 - t33 * t35 + t34 * t36, 0, 0, 0, 0, 0, t121 * t196 - t126 * t14 + t35 * t98 + t58 * t62, -t121 * t44 - t126 * t13 - t182 * t35 + t58 * t61, 0, 0, 0, 0, 0 (-qJD(6) * t201 - t13 * t168 + t172 * t35) * t97 + t200 * t62 + t14 * t63 - t196 * t26 -(qJD(6) * t200 + t13 * t172 + t168 * t35) * t97 - t201 * t62 + t14 * t65 - t196 * t25; 0, 0, 0, 0, t217 * t301, -t260 * t217, t215 * t229, -t215 * t230, 0, -t123 * t216 + (-pkin(9) * t216 * t254 + t166 * t181) * t174 + (-t166 * t209 + t180 * qJD(3) + (-qJD(3) * t166 + (-t166 ^ 2 - t159) * qJD(2)) * qJD(3) * pkin(2)) * t170 -(t170 * t181 + t236) * t166 + (-t159 * t289 + t278) * t252 + (pkin(9) * t230 - t308) * t216, -t121 * t82 - t122 * t81 + t290 * t130 - t131 * t34 - t132 * t33 + t291 * t133 - t22 * t135 + t136 * t21, -t21 * t81 + t22 * t82 + t290 * t34 - t291 * t33 + (t127 * t235 + (pkin(3) * t253 - t233) * t104) * t163, t108 * t61 - t182 * t79, -t107 * t61 - t108 * t62 + t182 * t80 - t79 * t98, t108 * t121 + t126 * t79 - t131 * t182 + t135 * t61, -t107 * t121 - t126 * t80 - t131 * t98 - t135 * t62, t121 * t135 + t126 * t131, t197 * t121 - t222 * t135 + t15 * t131 + t68 * t62 + t21 * t107 + t31 * t80 + t291 * t98 + t302 * t126 + (-t126 * t198 - t135 * t16) * qJD(5), t21 * t108 - t198 * t121 + t303 * t126 - t16 * t131 - t187 * t135 - t182 * t291 + t31 * t79 + t68 * t61, t25 * t89 + t27 * t65, -t25 * t88 - t26 * t89 - t27 * t63 - t28 * t65, t107 * t25 + t27 * t97 + t62 * t89 + t65 * t80, -t107 * t26 - t28 * t97 - t62 * t88 - t63 * t80, t107 * t62 + t80 * t97 (-t168 * t30 + t172 * t37) * t62 + t2 * t107 - t202 * t80 + t29 * t26 + t4 * t88 + t11 * t28 + (t168 * t211 - t172 * t206) * t97 + t293 * t63 -(t168 * t37 + t172 * t30) * t62 - t1 * t107 - t6 * t80 + t29 * t25 + t4 * t89 + t11 * t27 + (t168 * t206 + t172 * t211) * t97 + t293 * t65; 0, 0, 0, 0, -t275 * t301, t260 * t275, -t174 * t227, t170 * t227, 0 (-t203 * t269 + (-t180 - t232) * t170) * qJD(2) (-t125 * t273 + (t170 * t203 + t261) * t166) * qJD(2) + t261 * qJD(3) - t236 (t34 - t39) * t133 + (t33 - t40) * t130 + (-t121 * t162 - t122 * t165) * pkin(3), t33 * t39 - t34 * t40 + (-t104 * t231 + t162 * t22 - t165 * t21) * pkin(3), t169 * t61 - t182 * t218 (t61 - t288) * t173 + (-t62 + t280) * t169, t126 * t218 + t268 + t279, t194 + t286, -t126 * t133, -t15 * t133 + t158 * t62 - t39 * t98 + (-t21 + (-qJD(5) * t157 - t86) * t126) * t173 + (t40 * t126 + t184) * t169, t39 * t182 + t16 * t133 + t158 * t61 + t21 * t169 + (t157 * t248 + t292) * t126 + t184 * t173, t169 * t172 * t25 + (-t169 * t245 + t172 * t246 - t94) * t65, t63 * t94 + t65 * t93 + (-t168 * t65 - t172 * t63) * t246 + (-t284 - t172 * t26 + (t168 * t63 - t172 * t65) * qJD(6)) * t169, -t295 + t56 + (-t25 + t239) * t173 + (-t192 - t287) * t169, t296 + (t26 - t240) * t173 + (t193 - t304) * t169, -t173 * t62 + t277 * t97, t140 * t282 - t11 * t93 - t18 * t63 + (t168 * t219 + t172 * t281) * t97 + (t11 * t249 - t2 + (qJD(5) * t63 + t193) * t157) * t173 + (t11 * t244 + t202 * t130 + t157 * t26 + t4 * t168 + (t157 * t168 * t97 - t202) * qJD(5)) * t169, -t140 * t283 - t11 * t94 - t18 * t65 + (-t168 * t281 + t172 * t219) * t97 + (t11 * t247 + t1 + (qJD(5) * t65 + t192) * t157) * t173 + (-t11 * t245 + t6 * t130 + t157 * t25 + t4 * t172 + (t157 * t220 - t6) * qJD(5)) * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130 ^ 2 - t133 ^ 2, -t130 * t34 + t133 * t33 + t127, 0, 0, 0, 0, 0, t194 - t286, -t126 ^ 2 * t173 - t268 + t279, 0, 0, 0, 0, 0, t296 + (-t26 - t240) * t173 + (t193 + t304) * t169, t295 + t56 + (-t25 - t239) * t173 + (t192 - t287) * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182 * t98, t182 ^ 2 - t98 ^ 2, t61 + t288, -t62 - t280, t121, t182 * t31 - t222 + (-qJD(5) + t126) * t16, t126 * t15 + t31 * t98 - t187, t220 * t65 + t284 (t25 - t298) * t172 + (-t26 - t297) * t168, t182 * t65 + t220 * t97 + t283, -t168 * t97 ^ 2 - t182 * t63 + t282, t97 * t182, -pkin(5) * t26 - t16 * t63 + t183 * t168 - t299 * t172 - t182 * t202, -pkin(5) * t25 - t16 * t65 + t299 * t168 + t183 * t172 - t182 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t63, -t63 ^ 2 + t65 ^ 2, t25 + t298, -t26 + t297, t62, -t11 * t65 + t6 * t97 + t2, t11 * t63 - t202 * t97 - t1;];
tauc_reg  = t5;
