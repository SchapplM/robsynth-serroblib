% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRP10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRP10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:58:05
% EndTime: 2019-05-06 01:58:24
% DurationCPUTime: 7.65s
% Computational Cost: add. (18412->383), mult. (36067->476), div. (0->0), fcn. (24124->8), ass. (0->273)
t334 = pkin(7) + pkin(1);
t240 = sin(qJ(3));
t244 = cos(qJ(3));
t239 = sin(qJ(4));
t243 = cos(qJ(4));
t295 = qJD(1) * t244;
t216 = t239 * qJD(3) + t243 * t295;
t231 = t244 * qJDD(1);
t291 = qJD(1) * qJD(3);
t284 = t240 * t291;
t219 = t231 - t284;
t279 = -t243 * qJDD(3) + t239 * t219;
t182 = -t216 * qJD(4) - t279;
t214 = -t243 * qJD(3) + t239 * t295;
t264 = -t239 * qJDD(3) - t243 * t219;
t183 = -t214 * qJD(4) - t264;
t238 = sin(qJ(5));
t242 = cos(qJ(5));
t188 = t242 * t214 + t238 * t216;
t119 = -t188 * qJD(5) + t238 * t182 + t242 * t183;
t229 = t240 * qJD(1) + qJD(4);
t227 = qJD(5) + t229;
t318 = t188 * t227;
t348 = t119 - t318;
t283 = t244 * t291;
t289 = t240 * qJDD(1);
t218 = -t283 - t289;
t213 = qJDD(4) - t218;
t210 = qJDD(5) + t213;
t190 = -t238 * t214 + t242 * t216;
t317 = t190 * t188;
t134 = -t317 - t210;
t301 = t242 * t134;
t187 = t190 ^ 2;
t335 = t227 ^ 2;
t345 = -t187 - t335;
t100 = -t238 * t345 + t301;
t309 = t238 * t134;
t98 = t242 * t345 + t309;
t64 = t243 * t100 - t239 * t98;
t40 = t240 * t64 - t244 * t348;
t62 = t239 * t100 + t243 * t98;
t398 = -qJ(2) * t62 + t334 * t40;
t397 = pkin(3) * t62;
t396 = pkin(8) * t62;
t394 = pkin(3) * t348 - pkin(8) * t64;
t336 = t188 ^ 2;
t167 = t336 - t335;
t106 = -t238 * t167 + t301;
t110 = -t242 * t167 - t309;
t171 = t227 * t190;
t280 = t242 * t182 - t238 * t183;
t257 = t190 * qJD(5) - t280;
t89 = -t171 + t257;
t393 = -t240 * t89 + t244 * (t239 * t106 - t243 * t110);
t392 = pkin(4) * t98;
t391 = pkin(9) * t98;
t390 = pkin(9) * t100;
t344 = t187 - t336;
t346 = t171 + t257;
t54 = -t238 * t346 + t242 * t348;
t325 = t238 * t348;
t56 = t242 * t346 + t325;
t389 = -t240 * t344 + t244 * (t239 * t54 + t243 * t56);
t388 = t243 * t106 + t239 * t110;
t386 = t239 * t56 - t243 * t54;
t342 = -t317 + t210;
t308 = t238 * t342;
t341 = -t335 - t336;
t350 = t242 * t341 - t308;
t127 = t242 * t342;
t351 = t238 * t341 + t127;
t363 = t239 * t350 + t243 * t351;
t364 = -t239 * t351 + t243 * t350;
t377 = t240 * t364 - t244 * t346;
t384 = qJ(2) * t363 - t334 * t377;
t383 = pkin(3) * t363;
t382 = pkin(8) * t363;
t378 = -pkin(3) * t346 + pkin(8) * t364;
t168 = -t187 + t335;
t365 = t242 * t168 + t308;
t366 = -t238 * t168 + t127;
t376 = t239 * t366 + t243 * t365;
t347 = t119 + t318;
t375 = t244 * (-t239 * t365 + t243 * t366) + t240 * t347;
t121 = -t336 - t187;
t374 = pkin(3) * t121;
t373 = pkin(4) * t121;
t372 = pkin(4) * t351;
t371 = pkin(9) * t350;
t370 = pkin(9) * t351;
t367 = t244 * t121;
t195 = t216 * t214;
t343 = -t195 + t213;
t359 = t239 * t343;
t355 = t243 * t343;
t352 = t348 * qJ(6);
t247 = qJD(1) ^ 2;
t349 = t334 * t247;
t203 = t229 * t214;
t156 = t183 + t203;
t152 = (qJD(4) - t229) * t216 + t279;
t255 = (-t188 * t238 - t190 * t242) * t227;
t316 = t227 * t238;
t165 = t190 * t316;
t315 = t227 * t242;
t285 = t188 * t315;
t269 = t165 - t285;
t340 = t239 * t269 + t243 * t255;
t259 = t238 * t257 + t285;
t270 = t188 * t316 - t242 * t257;
t339 = t239 * t259 + t243 * t270;
t338 = t244 * (-t239 * t255 + t243 * t269) + t240 * t210;
t288 = t240 * t317;
t337 = t244 * (-t239 * t270 + t243 * t259) - t288;
t211 = t214 ^ 2;
t212 = t216 ^ 2;
t228 = t229 ^ 2;
t290 = qJD(2) * qJD(1);
t233 = 0.2e1 * t290;
t235 = qJDD(1) * qJ(2);
t241 = sin(qJ(1));
t245 = cos(qJ(1));
t271 = t245 * g(1) + t241 * g(2);
t261 = -t235 + t271;
t256 = t233 - t261;
t266 = -t219 + t284;
t267 = -t218 + t283;
t151 = t267 * pkin(3) + t266 * pkin(8) + t256 - t349;
t281 = t241 * g(1) - t245 * g(2);
t268 = qJDD(2) - t281;
t297 = t247 * qJ(2);
t254 = t268 - t297;
t202 = -t334 * qJDD(1) + t254;
t192 = t244 * g(3) - t240 * t202;
t246 = qJD(3) ^ 2;
t273 = pkin(3) * t240 - pkin(8) * t244;
t258 = t247 * t273;
t163 = -t246 * pkin(3) + qJDD(3) * pkin(8) - t240 * t258 - t192;
t113 = -t243 * t151 + t239 * t163;
t71 = pkin(4) * t343 - pkin(9) * t156 - t113;
t114 = t239 * t151 + t243 * t163;
t272 = t229 * pkin(4) - t216 * pkin(9);
t73 = -t211 * pkin(4) + t182 * pkin(9) - t229 * t272 + t114;
t44 = t238 * t73 - t242 * t71;
t45 = t238 * t71 + t242 * t73;
t20 = t238 * t45 - t242 * t44;
t333 = pkin(4) * t20;
t85 = t242 * t347;
t91 = (-qJD(5) + t227) * t190 + t280;
t55 = t238 * t91 - t85;
t332 = pkin(4) * t55;
t331 = pkin(5) * t242;
t330 = t257 * pkin(5);
t294 = qJD(6) * t227;
t221 = 0.2e1 * t294;
t142 = t188 * pkin(5) - t190 * qJ(6);
t278 = t210 * qJ(6) - t188 * t142 + t45;
t260 = -pkin(5) * t335 + t278;
t33 = t221 + t260;
t35 = -t210 * pkin(5) - qJ(6) * t335 + t190 * t142 + qJDD(6) + t44;
t329 = -pkin(5) * t35 + qJ(6) * t33;
t328 = -pkin(5) * t347 - qJ(6) * t89;
t326 = t238 * t347;
t324 = t239 * t20;
t321 = t243 * t20;
t320 = qJ(6) * t242;
t319 = qJDD(1) * pkin(1);
t314 = t229 * t239;
t313 = t229 * t243;
t236 = t240 ^ 2;
t312 = t236 * t247;
t237 = t244 ^ 2;
t311 = t237 * t247;
t191 = t240 * g(3) + t244 * t202;
t162 = qJDD(3) * pkin(3) + t246 * pkin(8) - t244 * t258 + t191;
t102 = t182 * pkin(4) + t211 * pkin(9) - t216 * t272 + t162;
t310 = t238 * t102;
t306 = t239 * t162;
t177 = t195 + t213;
t305 = t239 * t177;
t286 = t240 * t247 * t244;
t303 = t240 * (qJDD(3) + t286);
t302 = t242 * t102;
t300 = t243 * t162;
t299 = t243 * t177;
t298 = t244 * (qJDD(3) - t286);
t296 = t236 + t237;
t292 = qJD(4) + t229;
t287 = t240 * t195;
t282 = -qJ(6) * t238 - pkin(4);
t21 = t238 * t44 + t242 * t45;
t68 = t239 * t113 + t243 * t114;
t13 = t238 * t33 - t242 * t35;
t277 = pkin(4) * t13 + t329;
t53 = -t238 * t89 - t85;
t276 = pkin(4) * t53 + t328;
t275 = -t45 + t392;
t83 = t238 * t119 + t190 * t315;
t84 = t242 * t119 - t165;
t274 = t244 * (-t239 * t83 + t243 * t84) + t288;
t265 = t243 * t113 - t239 * t114;
t150 = t244 * t191 - t240 * t192;
t263 = -t44 + t372;
t262 = qJ(2) + t273;
t253 = -pkin(5) * t345 - qJ(6) * t134 + t260;
t252 = t253 - t392;
t251 = pkin(5) * t342 + qJ(6) * t341 - t35;
t250 = t251 + t372;
t249 = -pkin(5) * t171 + 0.2e1 * qJD(6) * t190 + t102;
t248 = t249 + t352;
t222 = t296 * qJDD(1);
t220 = t231 - 0.2e1 * t284;
t217 = 0.2e1 * t283 + t289;
t204 = -t254 + t319;
t201 = -t212 + t228;
t200 = t211 - t228;
t199 = t261 - 0.2e1 * t290 + t349;
t197 = -t303 + t244 * (-t246 - t311);
t196 = t240 * (-t246 - t312) + t298;
t194 = t212 - t211;
t193 = -t212 - t228;
t184 = -t228 - t211;
t175 = t211 + t212;
t157 = t292 * t214 + t264;
t155 = t183 - t203;
t153 = -t292 * t216 - t279;
t141 = -t239 * t193 - t299;
t140 = t243 * t193 - t305;
t136 = t243 * t184 - t359;
t135 = t239 * t184 + t355;
t116 = -t152 * t243 + t239 * t156;
t103 = t240 * t141 + t244 * t157;
t97 = t240 * t136 + t244 * t153;
t74 = t240 * t116 + t244 * t175;
t66 = -t302 - t391;
t61 = t244 * t162 + t240 * t68;
t60 = -t310 - t370;
t59 = t242 * t91 + t326;
t57 = -t242 * t89 + t326;
t51 = t239 * t84 + t243 * t83;
t42 = -pkin(4) * t348 - t310 + t390;
t39 = -pkin(4) * t346 + t302 + t371;
t38 = t248 - t330;
t31 = -t239 * t55 + t243 * t59;
t30 = -t239 * t53 + t243 * t57;
t29 = t239 * t59 + t243 * t55;
t28 = t239 * t57 + t243 * t53;
t27 = (-t257 - t346) * pkin(5) + t248;
t26 = t249 - t330 + 0.2e1 * t352;
t25 = -qJ(6) * t121 + t35;
t24 = t221 + (-t121 - t335) * pkin(5) + t278;
t23 = t240 * t31 - t367;
t22 = t240 * t30 - t367;
t19 = -t238 * t27 - t320 * t346 - t370;
t18 = -pkin(5) * t325 + t242 * t26 + t391;
t17 = pkin(4) * t102 + pkin(9) * t21;
t16 = t242 * t27 + t282 * t346 + t371;
t15 = -t390 + t238 * t26 + (pkin(4) + t331) * t348;
t14 = t238 * t35 + t242 * t33;
t12 = -pkin(9) * t55 - t20;
t11 = pkin(9) * t59 + t21 - t373;
t10 = -pkin(9) * t53 - t238 * t24 + t242 * t25;
t9 = t243 * t21 - t324;
t8 = t239 * t21 + t321;
t7 = pkin(9) * t57 + t238 * t25 + t242 * t24 - t373;
t6 = t244 * t102 + t240 * t9;
t5 = -pkin(9) * t13 + (-pkin(5) * t238 + t320) * t38;
t4 = -t239 * t13 + t243 * t14;
t3 = t243 * t13 + t239 * t14;
t2 = pkin(9) * t14 + (-t282 + t331) * t38;
t1 = t240 * t4 + t244 * t38;
t32 = [0, 0, 0, 0, 0, qJDD(1), t281, t271, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t268 - 0.2e1 * t319, t233 + 0.2e1 * t235 - t271, pkin(1) * t204 + qJ(2) * (-t247 * pkin(1) + t256), -t266 * t244, -t244 * t217 - t240 * t220, t298 - t240 * (t246 - t311), t267 * t240, t244 * (-t246 + t312) - t303, 0, qJ(2) * t217 - t334 * t196 - t240 * t199, qJ(2) * t220 - t334 * t197 - t244 * t199, t334 * t222 - t296 * t297 - t150, -qJ(2) * t199 - t334 * t150, t244 * (t243 * t183 - t216 * t314) + t287, t244 * (t243 * t153 - t239 * t155) + t240 * t194, t244 * (-t239 * t201 + t355) + t240 * t156, t244 * (-t239 * t182 + t214 * t313) - t287, t244 * (t243 * t200 - t305) - t240 * t152, t240 * t213 + t244 * (-t214 * t243 + t216 * t239) * t229, t244 * (-pkin(8) * t135 - t306) - t240 * (-pkin(3) * t135 + t113) + qJ(2) * t135 - t334 * t97, t244 * (-pkin(8) * t140 - t300) - t240 * (-pkin(3) * t140 + t114) + qJ(2) * t140 - t334 * t103, t244 * t265 - t334 * t74 + t262 * (-t152 * t239 - t243 * t156), -t262 * t265 - t334 * t61, t274, -t389, t375, t337, t393, t338, t244 * (-t239 * t39 + t243 * t60 - t382) - t240 * (-t263 - t383) + t384, t244 * (-t239 * t42 + t243 * t66 - t396) - t240 * (-t275 - t397) - t398, t244 * (-pkin(8) * t29 - t239 * t11 + t243 * t12) - t240 * (-pkin(3) * t29 - t332) + qJ(2) * t29 - t334 * t23, t244 * (-pkin(8) * t8 - pkin(9) * t321 - t239 * t17) - t240 * (-pkin(3) * t8 - t333) + qJ(2) * t8 - t334 * t6, t274, t375, t389, t338, -t393, t337, t244 * (-t239 * t16 + t243 * t19 - t382) - t240 * (-t250 - t383) + t384, t244 * (-pkin(8) * t28 + t243 * t10 - t239 * t7) - t240 * (-pkin(3) * t28 - t276) + qJ(2) * t28 - t334 * t22, t244 * (-t239 * t15 + t243 * t18 + t396) - t240 * (-t252 - 0.2e1 * t294 + t397) + t398, t244 * (-pkin(8) * t3 - t239 * t2 + t243 * t5) - t240 * (-pkin(3) * t3 - t277) + qJ(2) * t3 - t334 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t247, -t204, 0, 0, 0, 0, 0, 0, t196, t197, -t222, t150, 0, 0, 0, 0, 0, 0, t97, t103, t74, t61, 0, 0, 0, 0, 0, 0, t377, t40, t23, t6, 0, 0, 0, 0, 0, 0, t377, t22, -t40, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t286, (-t236 + t237) * t247, t231, -t286, -t289, qJDD(3), t191, t192, 0, 0, t239 * t183 + t216 * t313, t239 * t153 + t243 * t155, t243 * t201 + t359, t243 * t182 + t214 * t314, t239 * t200 + t299, (-t214 * t239 - t216 * t243) * t229, pkin(3) * t153 + pkin(8) * t136 + t300, pkin(3) * t157 + pkin(8) * t141 - t306, pkin(3) * t175 + pkin(8) * t116 + t68, pkin(3) * t162 + pkin(8) * t68, t51, -t386, t376, t339, -t388, t340, t239 * t60 + t243 * t39 + t378, t239 * t66 + t243 * t42 - t394, pkin(8) * t31 + t243 * t11 + t239 * t12 - t374, pkin(3) * t102 + pkin(8) * t9 - pkin(9) * t324 + t243 * t17, t51, t376, t386, t340, t388, t339, t243 * t16 + t239 * t19 + t378, pkin(8) * t30 + t239 * t10 + t243 * t7 - t374, t243 * t15 + t239 * t18 + t394, pkin(3) * t38 + pkin(8) * t4 + t243 * t2 + t239 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t194, t156, -t195, -t152, t213, -t113, -t114, 0, 0, t317, t344, t347, -t317, -t89, t210, t263, t275, t332, t333, t317, t347, -t344, t210, t89, -t317, t250, t276, t221 + t252, t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, t344, t347, -t317, -t89, t210, -t44, -t45, 0, 0, t317, t347, -t344, t210, t89, -t317, t251, t328, t221 + t253, t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342, t347, t345, t35;];
tauJ_reg  = t32;
