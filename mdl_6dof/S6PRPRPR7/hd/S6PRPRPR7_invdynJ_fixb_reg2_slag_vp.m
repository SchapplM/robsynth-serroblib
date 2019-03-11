% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:57
% EndTime: 2019-03-08 19:54:02
% DurationCPUTime: 3.68s
% Computational Cost: add. (3317->455), mult. (6687->572), div. (0->0), fcn. (4685->10), ass. (0->241)
t150 = -pkin(2) - pkin(8);
t141 = sin(pkin(6));
t145 = sin(qJ(2));
t273 = t141 * t145;
t235 = qJD(2) * t273;
t105 = qJD(1) * t235;
t148 = cos(qJ(2));
t247 = qJDD(1) * t141;
t227 = t148 * t247;
t178 = qJDD(3) + t105 - t227;
t142 = cos(pkin(6));
t262 = qJD(1) * t142;
t327 = qJD(4) * t262 - t150 * qJDD(2) - t178;
t147 = cos(qJ(4));
t144 = sin(qJ(4));
t261 = qJD(1) * t144;
t263 = qJD(1) * t141;
t232 = t148 * t263;
t198 = qJD(3) - t232;
t75 = t150 * qJD(2) + t198;
t44 = t142 * t261 - t147 * t75;
t326 = qJD(5) + t44;
t138 = t144 ^ 2;
t139 = t147 ^ 2;
t325 = t138 + t139;
t248 = qJD(2) * qJD(4);
t228 = t147 * t248;
t241 = t144 * qJDD(2);
t324 = t228 + t241;
t152 = qJD(2) ^ 2;
t268 = t148 * t152;
t181 = qJDD(2) * t145 + t268;
t272 = t141 * t148;
t240 = t144 * t272;
t252 = qJD(4) * t147;
t37 = -qJD(4) * t240 + t142 * t252 - t147 * t235;
t73 = t142 * t144 + t147 * t272;
t323 = -qJD(4) * t37 - qJDD(4) * t73 + (t144 * t181 + t145 * t228) * t141;
t134 = t147 * qJDD(2);
t229 = t144 * t248;
t175 = t229 - t134;
t36 = qJD(4) * t73 - t144 * t235;
t74 = t142 * t147 - t240;
t322 = -qJD(4) * t36 + qJDD(4) * t74 + (t145 * t175 - t147 * t268) * t141;
t174 = qJD(3) + (qJD(4) * pkin(9) - qJD(5)) * t147;
t284 = qJ(5) * t147;
t199 = pkin(9) * t144 - t284;
t180 = qJ(3) + t199;
t113 = t145 * t247;
t209 = t324 * pkin(4) + qJ(5) * t229 + t113;
t14 = t180 * qJDD(2) + (t174 + t232) * qJD(2) + t209;
t143 = sin(qJ(6));
t146 = cos(qJ(6));
t149 = -pkin(4) - pkin(9);
t258 = qJD(2) * t147;
t266 = pkin(5) * t258 + t326;
t20 = t149 * qJD(4) + t266;
t237 = t145 * t263;
t259 = qJD(2) * t144;
t201 = pkin(4) * t259 + t237;
t38 = qJD(2) * t180 + t201;
t195 = t143 * t38 - t146 * t20;
t246 = qJDD(1) * t142;
t254 = qJD(4) * t144;
t214 = t144 * t246 + t327 * t147 + t75 * t254;
t189 = qJDD(5) + t214;
t4 = -t175 * pkin(5) + t149 * qJDD(4) + t189;
t1 = -t195 * qJD(6) + t146 * t14 + t143 * t4;
t125 = qJD(6) + t258;
t321 = t125 * t195 + t1;
t9 = t143 * t20 + t146 * t38;
t2 = -qJD(6) * t9 - t143 * t14 + t146 * t4;
t320 = t9 * t125 + t2;
t231 = t146 * t259;
t255 = qJD(4) * t143;
t85 = -t231 + t255;
t186 = t125 * t85;
t250 = qJD(6) * t143;
t21 = qJD(4) * t250 - qJD(6) * t231 - t146 * qJDD(4) - t324 * t143;
t319 = t21 - t186;
t28 = -qJD(4) * pkin(4) + t326;
t253 = qJD(4) * t146;
t87 = t143 * t259 + t253;
t280 = qJD(6) * t87;
t22 = qJDD(4) * t143 - t324 * t146 + t280;
t295 = t125 * t87;
t318 = -t22 + t295;
t83 = -qJDD(6) + t175;
t66 = t146 * t83;
t317 = -t125 * t250 - t66;
t282 = qJD(4) * t85;
t316 = t282 + t66;
t256 = qJD(4) * qJ(5);
t236 = t147 * t262;
t45 = t144 * t75 + t236;
t31 = -t45 - t256;
t315 = pkin(4) * t252 + qJ(5) * t254;
t285 = cos(pkin(10));
t219 = t285 * t145;
t140 = sin(pkin(10));
t274 = t140 * t148;
t70 = t142 * t219 + t274;
t218 = t285 * t148;
t275 = t140 * t145;
t72 = -t142 * t275 + t218;
t314 = -g(1) * t72 - g(2) * t70;
t216 = qJ(3) - t284;
t51 = qJD(2) * t216 + t201;
t313 = t51 * t258 + qJDD(5);
t25 = t236 + (-pkin(5) * qJD(2) + t75) * t144;
t23 = t25 + t256;
t312 = t125 * t23 - t149 * t83;
t242 = qJDD(4) * t150;
t260 = qJD(2) * qJ(3);
t89 = t237 + t260;
t309 = qJD(4) * (t237 - t89 - t260) - t242;
t136 = t144 * pkin(4);
t92 = t136 + t216;
t308 = qJD(4) * (-qJD(2) * t92 + t237 - t51) - t242;
t215 = t327 * t144 - t147 * t246 - t75 * t252;
t158 = (t144 * t44 + t147 * t45) * qJD(4) - t215 * t144 - t214 * t147;
t278 = qJDD(4) * pkin(4);
t13 = t189 - t278;
t244 = qJDD(4) * qJ(5);
t165 = qJD(4) * qJD(5) - t215 + t244;
t159 = (t144 * t28 - t147 * t31) * qJD(4) + t165 * t144 - t13 * t147;
t202 = -t143 * t195 - t146 * t9;
t307 = -qJD(6) * t202 + t1 * t143 + t2 * t146;
t69 = -t142 * t218 + t275;
t306 = pkin(8) * t69;
t71 = t142 * t274 + t219;
t305 = pkin(8) * t71;
t270 = t145 * t147;
t169 = t141 * (-t143 * t270 + t146 * t148);
t78 = t136 + t180;
t297 = pkin(5) - t150;
t94 = t297 * t147;
t26 = -t143 * t78 + t146 * t94;
t48 = t174 + t315;
t221 = qJD(4) * t297;
t81 = t144 * t221;
t302 = -qJD(1) * t169 + qJD(6) * t26 - t143 * t81 + t146 * t48;
t269 = t146 * t147;
t170 = t141 * (-t143 * t148 - t145 * t269);
t27 = t143 * t94 + t146 * t78;
t301 = -qJD(1) * t170 - qJD(6) * t27 - t143 * t48 - t146 * t81;
t299 = t87 * t85;
t296 = t142 ^ 2 * qJDD(1) - g(3);
t294 = t143 * t22;
t293 = t143 * t83;
t292 = t144 * t21;
t291 = t144 * t87;
t197 = qJD(3) + t232;
t245 = qJDD(2) * qJ(3);
t47 = qJD(2) * t197 + t113 + t245;
t290 = t145 * t47;
t289 = t146 * t21;
t288 = t146 * t85;
t287 = t148 * t89;
t283 = qJD(2) * t89;
t281 = qJD(4) * t87;
t279 = qJDD(2) * pkin(2);
t277 = t125 * t143;
t276 = t140 * t141;
t271 = t144 * t152;
t267 = t23 * qJD(4);
t265 = pkin(2) * t272 + qJ(3) * t273;
t88 = pkin(4) * t258 + qJ(5) * t259;
t264 = t138 - t139;
t257 = qJD(2) * t148;
t251 = qJD(5) * t147;
t249 = qJD(6) * t146;
t243 = qJDD(4) * t144;
t238 = pkin(8) * t272 + t265;
t234 = t141 * t257;
t63 = t69 * pkin(2);
t226 = qJ(3) * t70 - t63;
t64 = t71 * pkin(2);
t225 = qJ(3) * t72 - t64;
t32 = t144 * t276 - t71 * t147;
t33 = t144 * t71 + t147 * t276;
t224 = -t32 * pkin(4) + qJ(5) * t33;
t220 = t141 * t285;
t34 = t144 * t220 + t69 * t147;
t35 = -t69 * t144 + t147 * t220;
t223 = t34 * pkin(4) - qJ(5) * t35;
t222 = -t73 * pkin(4) + qJ(5) * t74;
t95 = t325 * qJDD(2);
t213 = qJD(6) * t147 + qJD(2);
t151 = qJD(4) ^ 2;
t97 = qJDD(4) * t147 - t151 * t144;
t212 = t70 * t136 - t306 - t63;
t211 = t72 * t136 - t305 - t64;
t210 = t273 * t136 + t238;
t208 = t144 * t228;
t206 = g(1) * t71 + g(2) * t69;
t203 = t143 * t9 - t146 * t195;
t191 = t145 * (-qJD(2) * pkin(2) + t198) + t287;
t188 = t47 * qJ(3) + t89 * qJD(3);
t187 = t281 - t293;
t183 = -g(3) * t272 + t206;
t40 = t143 * t73 + t146 * t273;
t39 = -t143 * t273 + t146 * t73;
t179 = -t125 * t249 + t293;
t177 = g(1) * t32 - g(2) * t34 + g(3) * t73;
t176 = -g(1) * t33 + g(2) * t35 - g(3) * t74;
t171 = g(3) * t273 - t314;
t167 = -t113 + t171;
t166 = t183 + t227;
t164 = -t150 * t151 - t171;
t163 = t176 - t215;
t162 = t177 - t214;
t161 = qJDD(3) - t166;
t5 = -t324 * pkin(5) + t165;
t160 = -qJD(6) * t125 * t149 + t176 + t5;
t157 = t325 * t105 - t150 * t95 + t183;
t156 = qJD(4) * t45 + t162;
t19 = t216 * qJDD(2) + (t197 - t251) * qJD(2) + t209;
t61 = qJD(3) - t251 + t315;
t155 = -qJDD(2) * t92 - t19 + (-t61 + t232) * qJD(2) - t164;
t154 = t198 * qJD(2) + t164 + t245 + t47;
t153 = (-t144 * t74 + t147 * t73) * qJDD(2) + (t144 * t36 + t147 * t37 + (-t144 * t73 - t147 * t74) * qJD(4)) * qJD(2);
t120 = t147 * t271;
t101 = t264 * t152;
t96 = t147 * t151 + t243;
t93 = t297 * t144;
t82 = t147 * t221;
t80 = qJDD(2) * t139 - 0.2e1 * t208;
t79 = qJDD(2) * t138 + 0.2e1 * t208;
t77 = t97 - t271;
t76 = t243 + (t151 + t152) * t147;
t68 = t181 * t141;
t67 = (-qJDD(2) * t148 + t145 * t152) * t141;
t62 = pkin(9) * t258 + t88;
t53 = t178 - t279;
t52 = -0.2e1 * t144 * t134 + 0.2e1 * t264 * t248;
t18 = t143 * t25 + t146 * t62;
t17 = -t143 * t62 + t146 * t25;
t11 = qJD(6) * t39 + t143 * t37 + t146 * t234;
t10 = -qJD(6) * t40 - t143 * t234 + t146 * t37;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, -t67, -t68, 0 (t145 ^ 2 + t148 ^ 2) * t141 ^ 2 * qJDD(1) + t296, 0, 0, 0, 0, 0, 0, 0, t67, t68 (qJD(2) * t191 - t148 * t53 + t290) * t141 + t296, 0, 0, 0, 0, 0, 0, t323, -t322, t153, -t215 * t74 + t214 * t73 - t36 * t45 + t37 * t44 - g(3) + (t89 * t257 + t290) * t141, 0, 0, 0, 0, 0, 0, t153, -t323, t322, t165 * t74 + t13 * t73 + t28 * t37 + t31 * t36 - g(3) + (t145 * t19 + t51 * t257) * t141, 0, 0, 0, 0, 0, 0, t10 * t125 + t22 * t74 - t36 * t85 - t39 * t83, -t11 * t125 - t21 * t74 - t36 * t87 + t40 * t83, -t10 * t87 - t11 * t85 + t21 * t39 - t22 * t40, t1 * t40 - t10 * t195 + t11 * t9 + t2 * t39 - t23 * t36 + t5 * t74 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t166, t167, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t161 - 0.2e1 * t279, 0.2e1 * qJD(2) * qJD(3) - t167 + 0.2e1 * t245, -t53 * pkin(2) - g(1) * t225 - g(2) * t226 - g(3) * t265 - t191 * t263 + t188, t80, t52, t97, t79, -t96, 0, t154 * t144 - t309 * t147, t309 * t144 + t154 * t147, t157 - t158, -g(1) * (t225 - t305) - g(2) * (t226 - t306) - g(3) * t238 + (-t287 + (-t144 * t45 + t147 * t44) * t145) * t263 + t158 * t150 + t188, 0, -t97, t96, t80, t52, t79, t157 - t159, t155 * t144 + t308 * t147, -t308 * t144 + t155 * t147, t19 * t92 + t51 * t61 - g(1) * (t216 * t72 + t211) - g(2) * (t216 * t70 + t212) - g(3) * (-qJ(5) * t141 * t270 + t210) + (-t148 * t51 + (t144 * t31 + t147 * t28) * t145) * t263 + t159 * t150, t249 * t291 + (t87 * t252 - t292) * t143 (-t143 * t85 + t146 * t87) * t252 + (-t294 - t289 + (-t143 * t87 - t288) * qJD(6)) * t144 (t125 * t255 - t21) * t147 + (-t179 - t281) * t144, -t252 * t288 + (-t146 * t22 + t85 * t250) * t144 (t125 * t253 - t22) * t147 + (t282 + t317) * t144, -t125 * t254 - t147 * t83, -t93 * t22 - t26 * t83 - t82 * t85 + t206 * t146 + (-t143 * t314 - t23 * t253 + t2) * t147 - g(3) * t169 + t301 * t125 + (qJD(4) * t195 - t5 * t146 + t23 * t250 - t237 * t85) * t144, t93 * t21 + t27 * t83 - t82 * t87 - t206 * t143 + (-t146 * t314 + t23 * t255 - t1) * t147 - g(3) * t170 - t302 * t125 + (t9 * qJD(4) + t5 * t143 + t23 * t249 - t237 * t87) * t144, t21 * t26 - t22 * t27 - t301 * t87 - t302 * t85 - t202 * t252 + (-qJD(6) * t203 + t1 * t146 - t143 * t2 - t171) * t144, t1 * t27 + t2 * t26 - t5 * t93 - t23 * t82 - g(1) * (-pkin(5) * t71 + t211) - g(2) * (-pkin(5) * t69 + t212) - g(3) * t210 + t302 * t9 - t301 * t195 + (-g(3) * pkin(5) * t148 + (-g(3) * t199 - t23 * t261) * t145) * t141 + t314 * t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t152, t105 + t161 - t279 - t283, 0, 0, 0, 0, 0, 0, t77, -t76, -t95, t158 - t183 - t283, 0, 0, 0, 0, 0, 0, -t95, -t77, t76, -qJD(2) * t51 + t159 - t183, 0, 0, 0, 0, 0, 0, t144 * t22 + t316 * t147 + (t143 * t213 + t144 * t253) * t125, -t292 + t187 * t147 + (-t143 * t254 + t146 * t213) * t125 (-t87 * t254 + qJD(2) * t85 + (qJD(6) * t85 - t21) * t147) * t146 + (-t85 * t254 - qJD(2) * t87 + (t22 - t280) * t147) * t143, t202 * qJD(2) + (qJD(4) * t203 + t5) * t144 + (t267 - t307) * t147 - t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t101, t134, -t120, -t241, qJDD(4), -t89 * t258 + t156, -qJD(4) * t44 + t89 * t259 - t163, 0, 0, qJDD(4), -t134, t241, t120, -t101, -t120 (-pkin(4) * t147 - qJ(5) * t144) * qJDD(2), t88 * t259 - t156 - 0.2e1 * t278 + t313, 0.2e1 * t244 + (0.2e1 * qJD(5) + t44) * qJD(4) + (-t144 * t51 + t147 * t88) * qJD(2) + t163, -t13 * pkin(4) - g(1) * t224 - g(2) * t223 - g(3) * t222 + t165 * qJ(5) - t28 * t45 - t326 * t31 - t51 * t88, -t277 * t87 - t289 (-t22 - t295) * t146 + (t21 + t186) * t143 (-t147 * t277 + t291) * qJD(2) + t317, t146 * t186 + t294 (-t125 * t269 - t144 * t85) * qJD(2) + t179, t125 * t259, qJ(5) * t22 - t125 * t17 + t160 * t143 + t312 * t146 - t195 * t259 + t266 * t85, -qJ(5) * t21 + t125 * t18 - t312 * t143 + t160 * t146 - t9 * t259 + t266 * t87, t17 * t87 + t18 * t85 + (-t9 * t258 + t149 * t21 - t2 + (-t149 * t85 - t9) * qJD(6)) * t146 + (-t195 * t258 - t149 * t22 - t1 + (t149 * t87 - t195) * qJD(6)) * t143 + t177, t5 * qJ(5) - t9 * t18 + t195 * t17 - g(1) * (-pkin(9) * t32 + t224) - g(2) * (pkin(9) * t34 + t223) - g(3) * (-pkin(9) * t73 + t222) + t266 * t23 + t307 * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, qJDD(4) - t120, -t139 * t152 - t151, qJD(4) * t31 - t162 - t278 + t313, 0, 0, 0, 0, 0, 0, -t125 * t277 - t316, -t125 ^ 2 * t146 - t187, t318 * t143 + t319 * t146, t321 * t143 + t320 * t146 - t177 - t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, -t85 ^ 2 + t87 ^ 2, -t319, -t299, t318, -t83, -t23 * t87 - g(1) * (-t143 * t72 + t146 * t32) - g(2) * (-t143 * t70 - t146 * t34) - g(3) * t39 + t320, t23 * t85 - g(1) * (-t143 * t32 - t146 * t72) - g(2) * (t143 * t34 - t146 * t70) + g(3) * t40 - t321, 0, 0;];
tau_reg  = t3;
