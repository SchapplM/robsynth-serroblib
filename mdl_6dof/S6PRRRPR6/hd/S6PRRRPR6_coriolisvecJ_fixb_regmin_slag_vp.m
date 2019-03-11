% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:46
% EndTime: 2019-03-08 23:36:57
% DurationCPUTime: 4.07s
% Computational Cost: add. (3126->433), mult. (7749->591), div. (0->0), fcn. (5564->10), ass. (0->226)
t160 = sin(qJ(3));
t164 = cos(qJ(3));
t198 = pkin(3) * t160 - pkin(9) * t164;
t119 = t198 * qJD(3);
t124 = -t164 * pkin(3) - t160 * pkin(9) - pkin(2);
t159 = sin(qJ(4));
t161 = sin(qJ(2));
t163 = cos(qJ(4));
t243 = qJD(4) * t163;
t156 = sin(pkin(6));
t250 = qJD(1) * t156;
t165 = cos(qJ(2));
t258 = t164 * t165;
t314 = -(t159 * t161 + t163 * t258) * t250 + t159 * t119 + t124 * t243;
t236 = t164 * qJD(2);
t313 = qJD(4) - t236;
t301 = qJD(4) - qJD(6);
t247 = qJD(2) * t160;
t213 = t159 * t247;
t237 = t163 * qJD(3);
t110 = t213 - t237;
t239 = t159 * qJD(3);
t112 = t163 * t247 + t239;
t158 = sin(qJ(6));
t162 = cos(qJ(6));
t185 = -t162 * t110 + t158 * t112;
t54 = t158 * t110 + t162 * t112;
t312 = -t185 ^ 2 + t54 ^ 2;
t132 = t313 * qJD(5);
t234 = qJD(2) * qJD(3);
t149 = t160 * t234;
t140 = qJ(5) * t149;
t244 = qJD(4) * t159;
t216 = t161 * t250;
t121 = qJD(2) * pkin(8) + t216;
t157 = cos(pkin(6));
t248 = qJD(2) * t156;
t221 = t165 * t248;
t246 = qJD(3) * t160;
t249 = qJD(1) * t164;
t47 = -t121 * t246 + (qJD(3) * t157 + t221) * t249;
t265 = t157 * t160;
t138 = qJD(1) * t265;
t88 = t164 * t121 + t138;
t76 = qJD(3) * pkin(9) + t88;
t86 = (t119 + t216) * qJD(2);
t215 = t165 * t250;
t89 = qJD(2) * t124 - t215;
t178 = -t159 * t86 - t163 * t47 - t89 * t243 + t244 * t76;
t8 = t132 + t140 - t178;
t217 = t160 * t243;
t174 = t164 * t239 + t217;
t233 = qJD(3) * qJD(4);
t83 = qJD(2) * t174 + t159 * t233;
t6 = t83 * pkin(10) + t8;
t209 = -t159 * t47 + t163 * t86 - t76 * t243 - t89 * t244;
t297 = pkin(4) + pkin(5);
t212 = t164 * t237;
t82 = -qJD(2) * t212 + qJD(4) * t213 - t163 * t233;
t7 = t82 * pkin(10) - t297 * t149 - t209;
t223 = -t158 * t6 + t162 * t7;
t87 = -t160 * t121 + t157 * t249;
t199 = qJD(3) * pkin(3) + t87;
t175 = t112 * qJ(5) + t199;
t23 = -t297 * t110 + t175;
t311 = t23 * t54 - t223;
t134 = t313 * qJ(5);
t33 = t159 * t89 + t163 * t76;
t25 = t134 + t33;
t205 = pkin(4) * t149;
t9 = -t205 - t209;
t309 = -t25 * t313 + t9;
t308 = t23 * t185;
t307 = t54 * t185;
t259 = t163 * t164;
t146 = pkin(8) * t259;
t306 = qJD(4) * t146 - t163 * t119 + t124 * t244 - (t159 * t258 - t161 * t163) * t250;
t305 = qJ(5) * t246 + t314;
t304 = t159 * qJD(5) + t88;
t32 = -t159 * t76 + t163 * t89;
t254 = qJD(5) - t32;
t218 = t160 * t244;
t303 = t212 - t218;
t235 = -qJD(6) + t313;
t302 = -qJD(6) - t235;
t240 = qJD(6) * t162;
t241 = qJD(6) * t158;
t13 = t110 * t240 - t112 * t241 + t158 * t83 - t162 * t82;
t300 = t185 * t235 - t13;
t263 = t159 * qJ(5);
t299 = -t297 * t163 - t263;
t298 = t112 ^ 2;
t296 = pkin(9) - pkin(10);
t295 = pkin(10) * t160;
t224 = -pkin(8) * t159 - pkin(4);
t231 = pkin(10) * t259;
t294 = pkin(10) * t218 + (-t231 + (-pkin(5) + t224) * t160) * qJD(3) + t306;
t261 = t160 * t163;
t293 = -(-pkin(8) * qJD(3) + pkin(10) * qJD(4)) * t261 - (-qJD(5) + (-pkin(8) * qJD(4) + pkin(10) * qJD(3)) * t159) * t164 - t305;
t292 = -t164 * qJD(5) + (-t160 * t237 - t164 * t244) * pkin(8) + t305;
t291 = t224 * t246 + t306;
t275 = qJ(5) * t163;
t192 = pkin(4) * t159 - t275;
t290 = -t313 * t192 + t304;
t114 = t158 * t159 + t162 * t163;
t177 = t114 * t164;
t289 = -qJD(2) * t177 + t301 * t114;
t220 = t159 * t236;
t264 = t158 * t163;
t288 = t158 * t243 + t159 * t240 - t163 * t241 - t236 * t264 + (t220 - t244) * t162;
t179 = -t297 * t159 + t275;
t287 = t313 * t179 + t304;
t286 = qJD(2) * pkin(2);
t200 = t160 * t215;
t245 = qJD(3) * t164;
t48 = qJD(2) * t200 + qJD(3) * t138 + t121 * t245;
t173 = -t82 * qJ(5) + t112 * qJD(5) - t48;
t11 = t83 * pkin(4) - t173;
t285 = t11 * t159;
t284 = t11 * t163;
t22 = t110 * pkin(10) + t33;
t16 = t134 + t22;
t283 = t158 * t16;
t34 = t110 * pkin(4) - t175;
t281 = t34 * t112;
t280 = t48 * t159;
t279 = t48 * t163;
t278 = t199 * t159;
t277 = t82 * t159;
t116 = t198 * qJD(2);
t276 = t159 * t116 + t163 * t87;
t274 = t110 * qJ(5);
t273 = t110 * t313;
t272 = t112 * t110;
t271 = t112 * t313;
t270 = t313 * t159;
t269 = t313 * t163;
t268 = t156 * t161;
t267 = t156 * t165;
t168 = qJD(2) ^ 2;
t266 = t156 * t168;
t262 = t159 * t164;
t260 = t162 * t159;
t167 = qJD(3) ^ 2;
t257 = t167 * t160;
t256 = t167 * t164;
t255 = t112 * pkin(10) - t254;
t252 = t159 * t124 + t146;
t154 = t160 ^ 2;
t251 = -t164 ^ 2 + t154;
t242 = qJD(5) * t163;
t15 = -t297 * t313 - t255;
t232 = t15 * t240 + t158 * t7 + t162 * t6;
t230 = pkin(9) * t270;
t229 = pkin(9) * t269;
t228 = pkin(9) * t246;
t227 = pkin(9) * t237;
t36 = qJ(5) * t247 + t276;
t130 = t296 * t163;
t226 = t161 * t266;
t222 = t161 * t248;
t219 = t313 * t244;
t210 = -t158 * t82 - t162 * t83;
t208 = t163 * t116 - t159 * t87;
t145 = pkin(8) * t262;
t207 = t163 * t124 - t145;
t206 = t235 ^ 2;
t204 = t160 * t221;
t203 = t110 * t215;
t202 = t112 * t215;
t201 = t164 * t221;
t80 = -t164 * qJ(5) + t252;
t129 = t296 * t159;
t197 = pkin(10) * t220 - qJD(6) * t129 + t296 * t244 + t36;
t196 = (-t297 * t160 - t231) * qJD(2) - t208 - t301 * t130;
t122 = -t215 - t286;
t194 = -t122 - t215;
t193 = t163 * pkin(4) + t263;
t2 = t158 * t15 + t162 * t16;
t153 = t164 * pkin(4);
t49 = t164 * pkin(5) + t145 + t153 + (-t124 - t295) * t163;
t55 = t159 * t295 + t80;
t191 = t158 * t49 + t162 * t55;
t101 = t164 * t268 + t265;
t64 = t101 * t159 + t163 * t267;
t65 = t101 * t163 - t159 * t267;
t190 = -t65 * t158 + t64 * t162;
t189 = t64 * t158 + t65 * t162;
t24 = -pkin(4) * t313 + t254;
t188 = -t159 * t25 + t163 * t24;
t187 = t162 * qJ(5) - t158 * t297;
t186 = t158 * qJ(5) + t162 * t297;
t184 = -t260 + t264;
t183 = qJD(2) * t154 + t164 * t313;
t182 = pkin(8) + t192;
t181 = -t16 * t241 + t232;
t180 = t313 * t33 + t209;
t100 = -t157 * t164 + t160 * t268;
t95 = t114 * t160;
t176 = -pkin(8) + t179;
t172 = t313 * t32 + t178;
t171 = qJD(3) * (-t194 - t286);
t62 = -qJD(3) * t100 + t201;
t18 = qJD(4) * t65 + t62 * t159 - t163 * t222;
t63 = qJD(3) * t101 + t204;
t170 = t100 * t83 + t63 * t110 - t149 * t64 - t18 * t313;
t14 = qJD(6) * t54 + t210;
t19 = -qJD(4) * t64 + t159 * t222 + t62 * t163;
t169 = t100 * t82 - t63 * t112 + t149 * t65 + t19 * t313;
t123 = -pkin(3) - t193;
t108 = pkin(3) - t299;
t94 = t158 * t261 - t160 * t260;
t93 = t182 * t160;
t81 = t153 - t207;
t79 = t176 * t160;
t59 = t112 * pkin(4) + t274;
t41 = -t297 * t112 - t274;
t40 = -t82 + t273;
t39 = (qJD(4) * t193 - t242) * t160 + t182 * t245;
t38 = -pkin(4) * t247 - t208;
t29 = qJD(6) * t95 + t303 * t158 - t174 * t162;
t28 = t301 * t160 * t184 + qJD(3) * t177;
t26 = (t299 * qJD(4) + t242) * t160 + t176 * t245;
t10 = -t297 * t83 + t173;
t1 = t162 * t15 - t283;
t3 = [0, 0, -t226, -t165 * t266, 0, 0, 0, 0, 0, -t164 * t226 + (-t63 - t204) * qJD(3), t160 * t226 + (-t62 - t201) * qJD(3), 0, 0, 0, 0, 0, t170, -t169, t170, -t19 * t110 + t18 * t112 - t64 * t82 - t65 * t83, t169, t11 * t100 + t24 * t18 + t25 * t19 + t34 * t63 + t9 * t64 + t8 * t65, 0, 0, 0, 0, 0 -(-qJD(6) * t189 - t19 * t158 + t18 * t162) * t235 - t190 * t149 - t63 * t185 - t100 * t14 (qJD(6) * t190 + t18 * t158 + t19 * t162) * t235 + t189 * t149 - t63 * t54 - t100 * t13; 0, 0, 0, 0, 0.2e1 * t164 * t149, -0.2e1 * t251 * t234, t256, -t257, 0, -pkin(8) * t256 + t160 * t171, pkin(8) * t257 + t164 * t171, t303 * t112 - t82 * t261 (-t110 * t163 - t112 * t159) * t245 + (t277 - t163 * t83 + (t110 * t159 - t112 * t163) * qJD(4)) * t160, -t313 * t218 + t82 * t164 + (t112 * t160 + t163 * t183) * qJD(3), -t313 * t217 + t83 * t164 + (-t110 * t160 - t159 * t183) * qJD(3) (t313 - t236) * t246, -t306 * t313 + ((pkin(8) * t110 - t278) * qJD(3) - t209) * t164 + (-t203 - t199 * t243 + pkin(8) * t83 + t280 + (pkin(8) * t270 + qJD(2) * t207 + t32) * qJD(3)) * t160, -t314 * t313 + (-t199 * t237 + (qJD(3) * t112 + t219) * pkin(8) - t178) * t164 + (-t202 + t199 * t244 - pkin(8) * t82 + t279 + (pkin(8) * t269 - t252 * qJD(2) - t33) * qJD(3)) * t160, t39 * t110 + t93 * t83 + (t239 * t34 + t9) * t164 - t291 * t313 + (-t203 + t34 * t243 + t285 + (-qJD(2) * t81 - t24) * qJD(3)) * t160, -t80 * t83 - t81 * t82 + t291 * t112 - t292 * t110 + t188 * t245 + (-t159 * t8 + t163 * t9 + (-t159 * t24 - t163 * t25) * qJD(4)) * t160, -t39 * t112 + t93 * t82 + (-t237 * t34 - t8) * t164 + t292 * t313 + (t202 + t34 * t244 - t284 + (qJD(2) * t80 + t25) * qJD(3)) * t160, t11 * t93 + t8 * t80 + t9 * t81 + (t39 - t200) * t34 + t292 * t25 + t291 * t24, t13 * t95 + t54 * t28, -t13 * t94 - t95 * t14 - t185 * t28 - t54 * t29, t13 * t164 - t28 * t235 + (-qJD(2) * t95 - t54) * t246, t29 * t235 - t14 * t164 + (qJD(2) * t94 + t185) * t246 (t235 - t236) * t246, t223 * t164 + t26 * t185 + t79 * t14 + t10 * t94 + t23 * t29 - (t293 * t158 + t294 * t162) * t235 + (-t164 * t2 + t191 * t235) * qJD(6) + (t185 * t215 + (-(-t158 * t55 + t162 * t49) * qJD(2) - t1) * qJD(3)) * t160, -t181 * t164 + t26 * t54 + t79 * t13 + t10 * t95 + t23 * t28 - ((-qJD(6) * t49 + t293) * t162 + (qJD(6) * t55 - t294) * t158) * t235 + (t54 * t215 + (qJD(2) * t191 + t2) * qJD(3)) * t160; 0, 0, 0, 0, -t160 * t168 * t164, t251 * t168, 0, 0, 0, t88 * qJD(3) - t122 * t247 - t48, t194 * t236, t112 * t269 - t277 (-t82 - t273) * t163 + (-t271 - t83) * t159, t313 * t243 + (-t313 * t259 + (-t112 + t239) * t160) * qJD(2), -t219 + (t313 * t262 + (t110 + t237) * t160) * qJD(2), -t313 * t247, -pkin(3) * t83 - t279 - t208 * t313 - t88 * t110 + (-t229 - t278) * qJD(4) + (-t32 * t160 + (t164 * t199 - t228) * t159) * qJD(2), pkin(3) * t82 + t280 + t276 * t313 - t88 * t112 + (-t163 * t199 + t230) * qJD(4) + (t199 * t259 + (t33 - t227) * t160) * qJD(2), -t284 + t123 * t83 + t38 * t313 - t290 * t110 + (t159 * t34 - t229) * qJD(4) + (t160 * t24 + (-t164 * t34 - t228) * t159) * qJD(2), t36 * t110 - t38 * t112 + (t8 + t313 * t24 + (qJD(4) * t112 - t83) * pkin(9)) * t163 + ((qJD(4) * t110 - t82) * pkin(9) + t309) * t159, -t285 + t123 * t82 - t36 * t313 + t290 * t112 + (-t163 * t34 - t230) * qJD(4) + (t34 * t259 + (-t25 + t227) * t160) * qJD(2), t11 * t123 - t24 * t38 - t25 * t36 - t290 * t34 + (qJD(4) * t188 + t9 * t159 + t8 * t163) * pkin(9), -t13 * t184 + t289 * t54, -t13 * t114 + t14 * t184 - t185 * t289 - t288 * t54, -t289 * t235 + (qJD(3) * t184 + t54) * t247, t288 * t235 + (qJD(3) * t114 - t185) * t247, -t235 * t247, t10 * t114 + t108 * t14 + t287 * t185 + t288 * t23 - (t158 * t197 - t162 * t196) * t235 + (-(t162 * t129 - t158 * t130) * qJD(3) + t1) * t247, -t10 * t184 + t108 * t13 + t287 * t54 + t289 * t23 - (t158 * t196 + t162 * t197) * t235 + ((t158 * t129 + t162 * t130) * qJD(3) - t2) * t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, -t110 ^ 2 + t298, t40, -t83 + t271, t149, t112 * t199 + t180, -t110 * t199 + t172, -t59 * t110 + t180 + 0.2e1 * t205 - t281, pkin(4) * t82 - t83 * qJ(5) + (t25 - t33) * t112 + (t24 - t254) * t110, -t34 * t110 + t59 * t112 + 0.2e1 * t132 + 0.2e1 * t140 - t172, -t9 * pkin(4) + t8 * qJ(5) - t24 * t33 + t25 * t254 - t34 * t59, -t307, -t312, t300, t235 * t54 + t14, t149, t186 * t149 - t41 * t185 - (t158 * t255 - t162 * t22) * t235 + (t187 * t235 + t2) * qJD(6) + t311, t187 * t149 - t41 * t54 - t308 - (t158 * t22 + t162 * t255) * t235 + (-t186 * t235 - t283) * qJD(6) + t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149 + t272, t40, -t313 ^ 2 - t298, t281 + t309, 0, 0, 0, 0, 0, -t112 * t185 - t149 * t162 - t158 * t206, -t112 * t54 + t149 * t158 - t162 * t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, t312, -t300, t302 * t54 - t210, -t149, t302 * t2 - t311, -t1 * t235 - t181 + t308;];
tauc_reg  = t3;
