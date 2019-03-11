% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:38
% EndTime: 2019-03-08 22:40:51
% DurationCPUTime: 4.55s
% Computational Cost: add. (3690->425), mult. (10100->626), div. (0->0), fcn. (8056->12), ass. (0->223)
t163 = cos(pkin(7));
t167 = sin(qJ(3));
t272 = t163 * t167;
t155 = pkin(2) * t272;
t161 = sin(pkin(7));
t171 = cos(qJ(3));
t274 = t161 * t171;
t297 = pkin(4) + pkin(9);
t168 = sin(qJ(2));
t267 = t168 * t171;
t172 = cos(qJ(2));
t268 = t167 * t172;
t189 = t163 * t267 + t268;
t162 = sin(pkin(6));
t261 = qJD(1) * t162;
t84 = t189 * t261;
t309 = (t297 * t274 + t155) * qJD(3) - t84;
t270 = t167 * t168;
t242 = t162 * t270;
t219 = t163 * t242;
t254 = qJD(3) * t171;
t232 = t163 * t254;
t237 = t172 * t261;
t308 = -pkin(2) * t232 - qJD(1) * t219 + t171 * t237;
t255 = qJD(3) * t167;
t234 = t161 * t255;
t146 = pkin(3) * t234;
t205 = pkin(10) * t167 - qJ(4) * t171;
t253 = qJD(4) * t167;
t177 = (qJD(3) * t205 - t253) * t161;
t238 = t168 * t261;
t218 = t161 * t238;
t307 = t218 - t146 - t177;
t306 = -t163 * qJD(4) + t308;
t259 = qJD(2) * t161;
t113 = pkin(9) * t259 + t238;
t290 = qJD(2) * pkin(2);
t129 = t237 + t290;
t164 = cos(pkin(6));
t260 = qJD(1) * t164;
t239 = t161 * t260;
t52 = t167 * t113 + (-t129 * t163 - t239) * t171;
t265 = qJD(4) + t52;
t257 = qJD(2) * t163;
t150 = qJD(3) + t257;
t166 = sin(qJ(5));
t170 = cos(qJ(5));
t235 = t171 * t259;
t91 = t150 * t166 + t170 * t235;
t90 = qJD(6) + t91;
t258 = qJD(2) * t162;
t230 = qJD(1) * t258;
t305 = qJD(3) * t239 + t172 * t230;
t293 = -t297 * t234 - t306;
t256 = qJD(2) * t167;
t236 = t161 * t256;
t140 = qJD(5) + t236;
t165 = sin(qJ(6));
t169 = cos(qJ(6));
t216 = t166 * t235;
t93 = t150 * t170 - t216;
t59 = -t169 * t140 + t165 * t93;
t304 = t140 * t59;
t249 = qJD(5) * t170;
t251 = qJD(5) * t166;
t275 = t161 * t167;
t152 = pkin(9) * t275;
t240 = -pkin(2) * t171 - pkin(3);
t69 = pkin(4) * t275 + t152 + (-pkin(10) + t240) * t163;
t173 = -pkin(3) - pkin(10);
t280 = qJ(4) * t167;
t192 = t173 * t171 - t280;
t80 = (-pkin(2) + t192) * t161;
t303 = -t309 * t166 + t307 * t170 - t69 * t249 + t251 * t80;
t302 = t167 * t171;
t301 = pkin(9) * t274 + t155;
t266 = pkin(4) * t236 + t265;
t246 = qJD(2) * qJD(3);
t229 = t161 * t246;
t142 = t171 * t229;
t33 = t173 * t150 + t266;
t149 = t163 * t260;
t49 = t149 + (qJD(2) * t192 - t129) * t161;
t16 = t166 * t33 + t170 * t49;
t103 = t129 * t272;
t214 = t168 * t230;
t196 = t171 * t214;
t30 = qJD(3) * t103 + t113 * t254 + t163 * t196 + t305 * t167;
t26 = pkin(4) * t142 + t30;
t211 = t167 * t229;
t264 = pkin(3) * t211 + t161 * t214;
t48 = qJD(2) * t177 + t264;
t178 = -qJD(5) * t16 - t166 * t48 + t170 * t26;
t4 = -pkin(5) * t142 - t178;
t300 = (pkin(5) * t93 + t90 * pkin(11)) * t90 + t4;
t61 = t140 * t165 + t169 * t93;
t63 = -qJD(5) * t91 + t166 * t211;
t23 = qJD(6) * t61 - t169 * t142 + t165 * t63;
t282 = t301 * qJD(3) - t84;
t299 = -t282 * t150 - t30 * t163;
t12 = pkin(11) * t140 + t16;
t143 = t150 * qJ(4);
t53 = t171 * t113 + t167 * t239 + t103;
t45 = pkin(4) * t235 + t53;
t36 = t143 + t45;
t17 = pkin(5) * t91 - pkin(11) * t93 + t36;
t204 = t12 * t165 - t169 * t17;
t187 = t166 * t26 + t170 * t48 + t33 * t249 - t251 * t49;
t3 = pkin(11) * t142 + t187;
t197 = t167 * t214;
t207 = t113 * t255 - t129 * t232 + t163 * t197 - t305 * t171;
t27 = -t150 * qJD(4) + t207;
t20 = -pkin(4) * t211 - t27;
t64 = -qJD(5) * t216 + t150 * t249 - t170 * t211;
t8 = pkin(5) * t64 - pkin(11) * t63 + t20;
t1 = -qJD(6) * t204 + t165 * t8 + t169 * t3;
t198 = t166 * t69 + t170 * t80;
t298 = -qJD(5) * t198 + t307 * t166 + t309 * t170;
t296 = t59 * t90;
t295 = t61 * t90;
t233 = t161 * t254;
t294 = -pkin(5) * t233 - t298;
t148 = pkin(3) * t236;
t81 = t205 * t259 + t148;
t292 = t166 * t45 + t170 * t81;
t291 = pkin(9) * t234 + t306;
t289 = t140 * t93;
t288 = t165 * t64;
t287 = t169 * t64;
t225 = t169 * t90;
t286 = t173 * t90;
t247 = qJD(6) * t169;
t248 = qJD(6) * t165;
t22 = t140 * t247 + t165 * t142 + t169 * t63 - t248 * t93;
t285 = t22 * t165;
t283 = t91 * t140;
t209 = pkin(5) * t170 + pkin(11) * t166;
t281 = qJD(5) * t209 - (-pkin(4) - t209) * t236 + t265;
t279 = t140 * t166;
t278 = t140 * t170;
t277 = t140 * t173;
t158 = t161 ^ 2;
t174 = qJD(2) ^ 2;
t276 = t158 * t174;
t273 = t162 * t174;
t271 = t163 * t172;
t269 = t167 * t169;
t159 = t167 ^ 2;
t262 = -t171 ^ 2 + t159;
t252 = qJD(5) * t165;
t250 = qJD(5) * t169;
t245 = t158 * t290;
t244 = t167 * t278;
t243 = t166 * t274;
t241 = t168 * t273;
t94 = -t163 * qJ(4) - t301;
t231 = t140 * t249;
t132 = pkin(5) * t166 - pkin(11) * t170 + qJ(4);
t224 = pkin(11) * t235 - qJD(6) * t132 + t292;
t223 = t150 + t257;
t222 = 0.2e1 * t158 * t246;
t221 = t158 * t241;
t220 = t276 * t302;
t79 = pkin(4) * t274 - t94;
t217 = t161 * t168 * t258;
t108 = t163 * t166 + t170 * t274;
t109 = t163 * t170 - t243;
t37 = pkin(5) * t108 - pkin(11) * t109 + t79;
t215 = -pkin(11) * t233 - qJD(6) * t37 + t303;
t35 = pkin(11) * t275 + t198;
t75 = -qJD(5) * t108 + t166 * t234;
t76 = -qJD(5) * t243 + t163 * t249 - t170 * t234;
t210 = -t76 * pkin(5) + t75 * pkin(11) + qJD(6) * t35 - t293;
t206 = -pkin(3) * t171 - t280;
t6 = t12 * t169 + t165 * t17;
t107 = -t161 * t162 * t172 + t163 * t164;
t70 = -t164 * t274 + (-t171 * t271 + t270) * t162;
t47 = t107 * t170 + t166 * t70;
t188 = t163 * t268 + t267;
t71 = t162 * t188 + t164 * t275;
t203 = t165 * t71 + t169 * t47;
t202 = -t165 * t47 + t169 * t71;
t15 = -t166 * t49 + t170 * t33;
t201 = -t166 * t81 + t170 * t45;
t199 = -t166 * t80 + t170 * t69;
t46 = t107 * t166 - t170 * t70;
t195 = -t150 * t235 + t142;
t194 = -t247 * t90 - t288;
t193 = t248 * t90 - t287;
t191 = t167 * t36 + t173 * t254;
t190 = -t109 * t165 + t161 * t269;
t78 = t109 * t169 + t165 * t275;
t185 = t140 * t90;
t184 = t140 * t61;
t182 = t150 * t53 - t30;
t180 = (-qJ(4) * t254 - t253) * t161;
t11 = -pkin(5) * t140 - t15;
t179 = -pkin(11) * t64 + (t11 + t15) * t90;
t2 = -qJD(6) * t6 - t165 * t3 + t169 * t8;
t38 = t164 * t234 + (qJD(2) * t189 + qJD(3) * t188) * t162;
t176 = t107 * t211 - t150 * t38 - t171 * t221;
t39 = -qJD(2) * t219 - qJD(3) * t242 + (t172 * t258 + (t161 * t164 + t162 * t271) * qJD(3)) * t171;
t175 = t107 * t142 - t150 * t39 + t167 * t221;
t126 = t170 * t142;
t101 = -qJ(4) * t235 + t148;
t96 = t163 * t240 + t152;
t95 = (-pkin(2) + t206) * t161;
t88 = -t161 * t129 + t149;
t87 = (t165 * t171 + t166 * t269) * t259;
t86 = t165 * t166 * t236 - t169 * t235;
t83 = t146 + t180;
t57 = t149 + (qJD(2) * t206 - t129) * t161;
t54 = qJD(2) * t180 + t264;
t50 = t57 * t236;
t42 = -t143 - t53;
t41 = -pkin(3) * t150 + t265;
t34 = -pkin(5) * t275 - t199;
t32 = qJD(6) * t78 + t165 * t75 - t169 * t233;
t31 = qJD(6) * t190 + t165 * t233 + t169 * t75;
t18 = -pkin(5) * t235 - t201;
t14 = -qJD(5) * t46 + t38 * t166 + t170 * t217;
t13 = qJD(5) * t47 + t166 * t217 - t38 * t170;
t5 = [0, 0, -t241, -t172 * t273, 0, 0, 0, 0, 0, t176, t175 (t167 * t38 + t171 * t39 + (-t167 * t71 + t171 * t70) * qJD(3)) * t259, -t176, -t175, t107 * t54 + t217 * t57 - t27 * t71 + t30 * t70 + t38 * t41 - t39 * t42, 0, 0, 0, 0, 0, -t13 * t140 - t142 * t46 + t39 * t91 + t64 * t71, -t14 * t140 - t142 * t47 + t39 * t93 + t63 * t71, 0, 0, 0, 0, 0 (-qJD(6) * t203 - t14 * t165 + t39 * t169) * t90 + t202 * t64 + t13 * t59 + t46 * t23 -(qJD(6) * t202 + t14 * t169 + t39 * t165) * t90 - t203 * t64 + t13 * t61 + t46 * t22; 0, 0, 0, 0, t222 * t302, -t262 * t222, t223 * t233, -t223 * t234, 0 (t161 * t88 - t245) * t255 + t299, t207 * t163 + t308 * t150 + (-t171 * t245 + (pkin(9) * t150 * t167 + t171 * t88) * t161) * qJD(3) (t167 * t30 - t171 * t27 + (t167 * t42 + t171 * t41) * qJD(3) + ((qJD(3) * t96 - t291) * t171 + (qJD(3) * t94 + t282) * t167) * qJD(2)) * t161, -t158 * t196 + (-t57 * t255 + t171 * t54 + (t171 * t83 - t255 * t95) * qJD(2)) * t161 - t299, t158 * t197 - t27 * t163 - t291 * t150 + (-t57 * t254 - t167 * t54 + (-t167 * t83 - t254 * t95) * qJD(2)) * t161, t27 * t94 + t30 * t96 + t54 * t95 + (t83 - t218) * t57 + t291 * t42 + t282 * t41, t109 * t63 + t75 * t93, -t108 * t63 - t109 * t64 - t75 * t91 - t76 * t93, t75 * t140 + (t167 * t63 + (qJD(2) * t109 + t93) * t254) * t161, -t76 * t140 + (-t167 * t64 + (-qJD(2) * t108 - t91) * t254) * t161 (t140 * t161 + t158 * t256) * t254, t20 * t108 + t36 * t76 + t79 * t64 + t293 * t91 + t298 * t140 + (t178 * t167 + (qJD(2) * t199 + t15) * t254) * t161, t20 * t109 + t36 * t75 + t79 * t63 + t293 * t93 + t303 * t140 + (-t187 * t167 + (-t198 * qJD(2) - t16) * t254) * t161, t22 * t78 + t31 * t61, t190 * t22 - t23 * t78 - t31 * t59 - t32 * t61, t108 * t22 + t31 * t90 + t61 * t76 + t64 * t78, -t108 * t23 + t190 * t64 - t32 * t90 - t59 * t76, t108 * t64 + t76 * t90 (-t165 * t35 + t169 * t37) * t64 + t2 * t108 - t204 * t76 + t34 * t23 - t4 * t190 + t11 * t32 + (t165 * t215 - t169 * t210) * t90 + t294 * t59 -(t165 * t37 + t169 * t35) * t64 - t1 * t108 - t6 * t76 + t34 * t22 + t4 * t78 + t11 * t31 + (t165 * t210 + t169 * t215) * t90 + t294 * t61; 0, 0, 0, 0, -t220, t262 * t276, t195 (-qJD(3) + t150) * t236, 0, -t236 * t88 + t182, -t52 * t150 - t235 * t88 + t207 ((-qJ(4) * qJD(3) - t42 - t53) * t167 + (-pkin(3) * qJD(3) + t265 - t41) * t171) * t259, -t101 * t235 - t182 + t50, t265 * t150 + (t101 * t167 + t171 * t57) * t259 - t27, -t30 * pkin(3) - t27 * qJ(4) - t57 * t101 - t265 * t42 - t41 * t53, t63 * t170 - t279 * t93 (-t64 - t289) * t170 + (-t63 + t283) * t166, -t140 * t251 + t126 + (-t167 * t279 - t171 * t93) * t259, -t231 + (-t244 + (-qJD(3) * t166 + t91) * t171) * t259, -t140 * t235, qJ(4) * t64 + t20 * t166 - t201 * t140 + t266 * t91 + (-t166 * t277 + t170 * t36) * qJD(5) + (-t15 * t171 + t170 * t191) * t259, qJ(4) * t63 + t20 * t170 + t292 * t140 + t266 * t93 + (-t166 * t36 - t170 * t277) * qJD(5) + (t16 * t171 - t166 * t191) * t259, t22 * t169 * t170 + (-t166 * t250 - t170 * t248 - t87) * t61, t87 * t59 + t61 * t86 + (t165 * t61 + t169 * t59) * t251 + (-t285 - t169 * t23 + (t165 * t59 - t169 * t61) * qJD(6)) * t170, -t87 * t90 + (-t250 * t90 + t22) * t166 + (t184 - t193) * t170, t86 * t90 + (t252 * t90 - t23) * t166 + (t194 - t304) * t170, t64 * t166 + t90 * t278, t132 * t287 - t11 * t86 - t18 * t59 + (t224 * t165 + t281 * t169) * t90 + (-t11 * t252 + t2 + (qJD(5) * t59 + t194) * t173) * t166 + (-t204 * t236 + t11 * t247 + t4 * t165 - t173 * t23 + (-t165 * t286 - t204) * qJD(5)) * t170, -t132 * t288 - t11 * t87 - t18 * t61 + (-t281 * t165 + t224 * t169) * t90 + (-t11 * t250 - t1 + (qJD(5) * t61 + t193) * t173) * t166 + (-t6 * t236 - t11 * t248 + t4 * t169 - t173 * t22 + (-t169 * t286 - t6) * qJD(5)) * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t220, -t150 ^ 2 - t159 * t276, t150 * t42 + t30 + t50, 0, 0, 0, 0, 0, -t140 * t279 - t150 * t91 + t126, -t231 - t150 * t93 + (-t166 * t254 - t244) * t259, 0, 0, 0, 0, 0, -t150 * t225 + (-t165 * t185 - t23) * t170 + (t194 + t304) * t166, t165 * t150 * t90 + (-t169 * t185 - t22) * t170 + (t184 + t193) * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93 * t91, -t91 ^ 2 + t93 ^ 2, t63 + t283, -t64 + t289, t142, t16 * t140 - t36 * t93 + t178, t140 * t15 + t36 * t91 - t187, t225 * t61 + t285 (t22 - t296) * t169 + (-t23 - t295) * t165, t225 * t90 - t61 * t93 + t288, -t90 ^ 2 * t165 + t59 * t93 + t287, -t90 * t93, -pkin(5) * t23 - t16 * t59 + t179 * t165 - t300 * t169 + t204 * t93, -pkin(5) * t22 - t16 * t61 + t300 * t165 + t179 * t169 + t6 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t59, -t59 ^ 2 + t61 ^ 2, t22 + t296, -t23 + t295, t64, -t11 * t61 + t6 * t90 + t2, t11 * t59 - t204 * t90 - t1;];
tauc_reg  = t5;
