% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:23
% EndTime: 2021-01-16 03:28:44
% DurationCPUTime: 6.24s
% Computational Cost: add. (5720->450), mult. (13725->626), div. (0->0), fcn. (11451->18), ass. (0->243)
t202 = cos(qJ(6));
t274 = qJD(6) * t202;
t192 = sin(pkin(12));
t195 = cos(pkin(12));
t204 = cos(qJ(3));
t285 = t195 * t204;
t260 = qJD(2) * t285;
t200 = sin(qJ(3));
t278 = qJD(2) * t200;
t144 = t192 * t278 - t260;
t203 = cos(qJ(5));
t131 = t203 * t144;
t153 = t192 * t204 + t195 * t200;
t147 = t153 * qJD(2);
t199 = sin(qJ(5));
t87 = -t147 * t199 - t131;
t342 = t202 * t87;
t350 = t274 - t342;
t188 = qJD(3) + qJD(5);
t306 = t188 * t87;
t146 = t153 * qJD(3);
t269 = t204 * qJDD(2);
t270 = t200 * qJDD(2);
t235 = t192 * t270 - t195 * t269;
t100 = qJD(2) * t146 + t235;
t272 = qJD(2) * qJD(3);
t258 = t200 * t272;
t218 = qJDD(2) * t153 - t192 * t258;
t257 = t204 * t272;
t101 = t195 * t257 + t218;
t276 = qJD(5) * t199;
t36 = -qJD(5) * t131 - t199 * t100 + t203 * t101 - t147 * t276;
t349 = t36 - t306;
t197 = qJ(4) + pkin(8);
t251 = qJD(3) * t197;
t137 = qJD(4) * t204 - t200 * t251;
t138 = -qJD(4) * t200 - t204 * t251;
t194 = sin(pkin(6));
t205 = cos(qJ(2));
t287 = t194 * t205;
t263 = qJD(1) * t287;
t297 = -t137 * t192 + t195 * t138 + t153 * t263;
t152 = t192 * t200 - t285;
t296 = t195 * t137 + t192 * t138 + t152 * t263;
t233 = t144 * t199 - t203 * t147;
t212 = qJD(5) * t233 - t203 * t100 - t101 * t199;
t307 = t188 * t233;
t348 = t212 - t307;
t284 = -qJD(6) + t87;
t347 = t284 + qJD(6);
t198 = sin(qJ(6));
t187 = qJDD(3) + qJDD(5);
t275 = qJD(6) * t198;
t24 = t198 * t187 + t188 * t274 + t202 * t36 + t233 * t275;
t21 = t24 * t202;
t73 = t188 * t198 - t202 * t233;
t25 = qJD(6) * t73 - t202 * t187 + t198 * t36;
t71 = -t202 * t188 - t198 * t233;
t346 = -t198 * t25 - t350 * t71 + t21;
t20 = t24 * t198;
t345 = t350 * t73 + t20;
t34 = qJDD(6) - t212;
t31 = t198 * t34;
t311 = -t274 * t284 + t31;
t315 = t73 * t233;
t344 = t284 * t342 + t311 + t315;
t320 = pkin(9) * t147;
t201 = sin(qJ(2));
t279 = qJD(1) * t201;
t264 = t194 * t279;
t241 = t197 * qJD(2) + t264;
t196 = cos(pkin(6));
t280 = qJD(1) * t196;
t112 = t200 * t280 + t204 * t241;
t104 = t192 * t112;
t111 = -t200 * t241 + t204 * t280;
t309 = qJD(3) * pkin(3);
t108 = t111 + t309;
t59 = t195 * t108 - t104;
t43 = qJD(3) * pkin(4) - t320 + t59;
t321 = pkin(9) * t144;
t286 = t195 * t112;
t60 = t192 * t108 + t286;
t46 = t60 - t321;
t22 = -t199 * t46 + t203 * t43;
t17 = -pkin(5) * t188 - t22;
t343 = t17 * t87;
t295 = cos(pkin(11));
t248 = t295 * t205;
t193 = sin(pkin(11));
t291 = t193 * t201;
t140 = t196 * t291 - t248;
t249 = t295 * t201;
t290 = t193 * t205;
t142 = t196 * t249 + t290;
t189 = qJ(3) + pkin(12);
t186 = qJ(5) + t189;
t178 = sin(t186);
t179 = cos(t186);
t250 = t194 * t295;
t289 = t194 * t201;
t292 = t193 * t194;
t223 = -g(3) * (-t178 * t289 + t179 * t196) - g(2) * (-t142 * t178 - t179 * t250) - g(1) * (t140 * t178 + t179 * t292);
t271 = t196 * qJDD(1);
t168 = t204 * t271;
t273 = qJD(1) * qJD(2);
t123 = qJDD(2) * pkin(8) + (qJDD(1) * t201 + t205 * t273) * t194;
t213 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t280 + t123;
t231 = t241 * qJD(3);
t54 = qJDD(3) * pkin(3) - t200 * t213 - t204 * t231 + t168;
t55 = (-t231 + t271) * t200 + t213 * t204;
t29 = -t192 * t55 + t195 * t54;
t13 = qJDD(3) * pkin(4) - pkin(9) * t101 + t29;
t30 = t192 * t54 + t195 * t55;
t14 = -pkin(9) * t100 + t30;
t23 = t199 * t43 + t203 * t46;
t330 = -qJD(5) * t23 + t203 * t13 - t199 * t14;
t3 = -pkin(5) * t187 - t330;
t222 = t223 - t3;
t341 = t233 * t87;
t316 = t71 * t233;
t149 = t152 * qJD(3);
t340 = -pkin(9) * t149 - t297;
t339 = -pkin(9) * t146 + t296;
t338 = t198 * t284;
t337 = t284 * t233;
t141 = -t196 * t248 + t291;
t143 = t196 * t290 + t249;
t238 = g(1) * t143 + g(2) * t141;
t220 = g(3) * t287 - t238;
t215 = t220 * t179;
t102 = t203 * t152 + t153 * t199;
t103 = -t152 * t199 + t153 * t203;
t181 = pkin(3) * t204 + pkin(2);
t125 = pkin(4) * t152 - t181;
t39 = pkin(5) * t102 - pkin(10) * t103 + t125;
t336 = t39 * t34 - t215;
t335 = t233 ^ 2 - t87 ^ 2;
t121 = t178 * t196 + t179 * t289;
t327 = (qJD(5) * t43 + t14) * t203 + t199 * t13 - t46 * t276;
t93 = t142 * t179 - t178 * t250;
t95 = -t140 * t179 + t178 * t292;
t134 = -t181 * qJD(2) + qJD(4) - t263;
t96 = pkin(4) * t144 + t134;
t334 = g(1) * t95 + g(2) * t93 + g(3) * t121 - t96 * t87 - t327;
t18 = pkin(10) * t188 + t23;
t35 = -pkin(5) * t87 + pkin(10) * t233 + t96;
t5 = t18 * t202 + t198 * t35;
t333 = t17 * t274 - t222 * t198 - t5 * t233;
t234 = t18 * t198 - t202 * t35;
t332 = t17 * t275 - t234 * t233;
t331 = t233 * t96 + t223 + t330;
t52 = -pkin(5) * t233 - pkin(10) * t87;
t32 = t202 * t34;
t227 = -t275 * t284 - t32;
t183 = t200 * t309;
t236 = pkin(4) * t146 + t183 - t264;
t323 = pkin(3) * t195;
t180 = pkin(4) + t323;
t324 = pkin(3) * t192;
t282 = t199 * t180 + t203 * t324;
t2 = pkin(10) * t187 + t327;
t239 = g(1) * t140 - g(2) * t142;
t221 = -g(3) * t289 + t239;
t160 = t197 * t200;
t161 = t197 * t204;
t116 = -t195 * t160 - t161 * t192;
t77 = -pkin(9) * t153 + t116;
t117 = -t192 * t160 + t195 * t161;
t78 = -pkin(9) * t152 + t117;
t41 = t199 * t78 - t203 * t77;
t314 = qJD(5) * t41 + t340 * t199 - t339 * t203;
t42 = t199 * t77 + t203 * t78;
t56 = -qJD(5) * t102 - t146 * t199 - t149 * t203;
t329 = -(qJD(6) * t35 + t2) * t102 + t17 * t56 + t3 * t103 - (-qJD(6) * t39 + t314) * t284 - t42 * t34 + t221;
t259 = t201 * t273;
t165 = t194 * t259;
t206 = qJD(3) ^ 2;
t255 = qJDD(1) * t287;
t328 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t206 + t194 * (-g(3) * t205 + t259) - t165 + t238 + t255;
t322 = pkin(3) * t200;
t150 = t196 * t204 - t200 * t289;
t318 = g(3) * t150;
t313 = qJD(5) * t42 + t339 * t199 + t340 * t203;
t310 = qJD(2) * pkin(2);
t308 = t103 * t17;
t305 = t198 * t73;
t303 = t202 * t73;
t302 = t202 * t284;
t229 = t180 * t203 - t199 * t324;
t62 = -t111 * t192 - t286;
t49 = t62 + t321;
t64 = t195 * t111 - t104;
t50 = t64 - t320;
t299 = -t229 * qJD(5) + t199 * t49 + t203 * t50;
t298 = t282 * qJD(5) - t199 * t50 + t203 * t49;
t288 = t194 * t204;
t283 = qJDD(1) - g(3);
t190 = t200 ^ 2;
t281 = -t204 ^ 2 + t190;
t277 = qJD(2) * t201;
t182 = pkin(3) * t278;
t266 = t198 * t287;
t265 = t202 * t287;
t262 = t194 * t277;
t261 = qJD(2) * t287;
t256 = t205 * t272;
t118 = pkin(4) * t147 + t182;
t244 = t194 * t283;
t136 = pkin(10) + t282;
t242 = qJD(6) * t136 + t118 + t52;
t57 = qJD(5) * t103 + t203 * t146 - t149 * t199;
t237 = pkin(5) * t57 - pkin(10) * t56 + t236;
t151 = t196 * t200 + t201 * t288;
t82 = t150 * t195 - t151 * t192;
t83 = t150 * t192 + t151 * t195;
t44 = t199 * t83 - t203 * t82;
t45 = t199 * t82 + t203 * t83;
t232 = -t338 * t87 - t227;
t207 = qJD(2) ^ 2;
t230 = qJDD(2) * t205 - t201 * t207;
t228 = -g(1) * t193 + t295 * g(2);
t226 = -t198 * t45 - t265;
t225 = -t202 * t45 + t266;
t159 = -t263 - t310;
t219 = -qJD(2) * t159 - t123 - t239;
t217 = -t136 * t34 - t284 * t299 - t343;
t216 = pkin(3) * t258 - t181 * qJDD(2) + qJDD(4) + t165;
t211 = -pkin(8) * qJDD(3) + (t159 + t263 - t310) * qJD(3);
t99 = t216 - t255;
t58 = pkin(4) * t100 + t99;
t185 = cos(t189);
t184 = sin(t189);
t135 = -pkin(5) - t229;
t110 = -qJD(3) * t151 - t200 * t261;
t109 = qJD(3) * t150 + t204 * t261;
t63 = t109 * t195 + t110 * t192;
t61 = -t109 * t192 + t110 * t195;
t9 = qJD(5) * t45 + t199 * t63 - t203 * t61;
t8 = -qJD(5) * t44 + t199 * t61 + t203 * t63;
t7 = -pkin(5) * t212 - pkin(10) * t36 + t58;
t6 = t202 * t7;
t1 = [t283, 0, t230 * t194, (-qJDD(2) * t201 - t205 * t207) * t194, 0, 0, 0, 0, 0, qJD(3) * t110 + qJDD(3) * t150 + (-t200 * t256 + t204 * t230) * t194, -qJD(3) * t109 - qJDD(3) * t151 + (-t200 * t230 - t204 * t256) * t194, qJD(3) * t61 + qJDD(3) * t82 + (-t100 * t205 + t144 * t277) * t194, -qJD(3) * t63 - qJDD(3) * t83 + (-t101 * t205 + t147 * t277) * t194, -t100 * t83 - t101 * t82 - t144 * t63 - t147 * t61, t29 * t82 + t30 * t83 + t59 * t61 + t60 * t63 - g(3) + (t134 * t277 - t205 * t99) * t194, 0, 0, 0, 0, 0, -t187 * t44 - t188 * t9 + (t205 * t212 - t277 * t87) * t194, -t187 * t45 - t188 * t8 + (-t205 * t36 - t233 * t277) * t194, 0, 0, 0, 0, 0, -(qJD(6) * t225 - t198 * t8 + t202 * t262) * t284 + t226 * t34 + t9 * t71 + t44 * t25, (qJD(6) * t226 + t198 * t262 + t202 * t8) * t284 + t225 * t34 + t9 * t73 + t44 * t24; 0, qJDD(2), t283 * t287 + t238, -t201 * t244 - t239, qJDD(2) * t190 + 0.2e1 * t200 * t257, 0.2e1 * t200 * t269 - 0.2e1 * t281 * t272, qJDD(3) * t200 + t204 * t206, qJDD(3) * t204 - t200 * t206, 0, t211 * t200 + t328 * t204, -t328 * t200 + t211 * t204, -t144 * t264 + qJDD(3) * t116 - t100 * t181 + t134 * t146 + t152 * t99 - t220 * t185 + (t144 * t322 + t297) * qJD(3), -t147 * t264 - qJDD(3) * t117 - t101 * t181 - t134 * t149 + t153 * t99 + t220 * t184 + (t147 * t322 - t296) * qJD(3), -t100 * t117 - t101 * t116 - t296 * t144 - t146 * t60 - t297 * t147 + t149 * t59 - t152 * t30 - t153 * t29 + t221, t30 * t117 + t29 * t116 - t99 * t181 + t134 * t183 - g(1) * (-t140 * t197 - t143 * t181) - g(2) * (-t141 * t181 + t142 * t197) + t296 * t60 + t297 * t59 + (-t134 * t279 - g(3) * (t181 * t205 + t197 * t201)) * t194, t103 * t36 - t233 * t56, -t102 * t36 + t103 * t212 + t233 * t57 + t56 * t87, t103 * t187 + t188 * t56, -t102 * t187 - t188 * t57, 0, t102 * t58 - t125 * t212 - t187 * t41 - t188 * t313 - t236 * t87 + t57 * t96 - t215, t103 * t58 + t125 * t36 + t178 * t220 - t187 * t42 + t188 * t314 - t233 * t236 + t56 * t96, t56 * t303 + (-t275 * t73 + t21) * t103, (-t202 * t71 - t305) * t56 + (-t20 - t202 * t25 + (t198 * t71 - t303) * qJD(6)) * t103, t102 * t24 - t103 * t227 - t302 * t56 + t57 * t73, -t102 * t25 - t103 * t311 + t338 * t56 - t57 * t71, t102 * t34 - t284 * t57, t6 * t102 + t41 * t25 - t234 * t57 + t313 * t71 + (-t237 * t284 + (-t102 * t18 + t284 * t42 + t308) * qJD(6) + t336) * t202 + t329 * t198, t41 * t24 - t5 * t57 + t313 * t73 + (-(-qJD(6) * t18 + t7) * t102 - qJD(6) * t308 - (qJD(6) * t42 - t237) * t284 - t336) * t198 + t329 * t202; 0, 0, 0, 0, -t200 * t207 * t204, t281 * t207, t270, t269, qJDD(3), t200 * t219 + t228 * t288 + t168 - t318, g(3) * t151 + (-t194 * t228 - t271) * t200 + t219 * t204, -t62 * qJD(3) - t134 * t147 - g(1) * (t140 * t184 + t185 * t292) - g(2) * (-t142 * t184 - t185 * t250) - g(3) * (-t184 * t289 + t185 * t196) + (qJDD(3) * t195 - t144 * t278) * pkin(3) + t29, t64 * qJD(3) + t134 * t144 - g(1) * (t140 * t185 - t184 * t292) - g(2) * (-t142 * t185 + t184 * t250) - g(3) * (-t184 * t196 - t185 * t289) + (-qJDD(3) * t192 - t147 * t278) * pkin(3) - t30, (t60 + t62) * t147 + (-t59 + t64) * t144 + (-t100 * t192 - t101 * t195) * pkin(3), -t134 * t182 + t29 * t323 + t30 * t324 - t59 * t62 - t60 * t64 + (-g(1) * (t140 * t200 + t193 * t288) - g(2) * (-t142 * t200 - t204 * t250) - t318) * pkin(3), t341, t335, t349, t348, t187, t118 * t87 + t187 * t229 - t188 * t298 + t331, t118 * t233 - t187 * t282 + t188 * t299 + t334, t345, t284 * t305 + t346, t344, t232 - t316, -t337, t135 * t25 + t298 * t71 + t217 * t198 + (t242 * t284 + t222) * t202 + t332, t135 * t24 + t202 * t217 - t242 * t338 + t298 * t73 + t333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t147 + t235, (-t144 + t260) * qJD(3) + t218, -t144 ^ 2 - t147 ^ 2, t144 * t60 + t147 * t59 - t205 * t244 + t216 - t238, 0, 0, 0, 0, 0, -t212 - t307, t36 + t306, 0, 0, 0, 0, 0, t232 + t316, -t284 * t302 - t31 + t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t341, t335, t349, t348, t187, t188 * t23 + t331, t188 * t22 + t334, t345, t338 * t73 + t346, t344, -t284 * t338 - t316 + t32, -t337, -pkin(5) * t25 - t23 * t71 + (-pkin(10) * t34 - t22 * t284 - t343) * t198 + (-(-pkin(10) * qJD(6) - t52) * t284 + t222) * t202 + t332, -pkin(5) * t24 - (t198 * t52 + t202 * t22) * t284 - t23 * t73 - t17 * t342 + t227 * pkin(10) + t333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, -t284 * t71 + t24, -t284 * t73 - t25, t34, -t198 * t2 + t6 - t17 * t73 - g(1) * (t143 * t202 - t198 * t95) - g(2) * (t141 * t202 - t198 * t93) - g(3) * (-t121 * t198 - t265) - t347 * t5, -t202 * t2 - t198 * t7 + t17 * t71 - g(1) * (-t143 * t198 - t202 * t95) - g(2) * (-t141 * t198 - t202 * t93) - g(3) * (-t121 * t202 + t266) + t347 * t234;];
tau_reg = t1;
