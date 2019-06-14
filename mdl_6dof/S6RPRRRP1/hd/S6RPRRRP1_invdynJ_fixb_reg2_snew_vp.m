% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:09:19
% EndTime: 2019-05-06 01:09:38
% DurationCPUTime: 7.23s
% Computational Cost: add. (19311->435), mult. (38181->570), div. (0->0), fcn. (26335->10), ass. (0->270)
t259 = sin(qJ(4));
t262 = cos(qJ(4));
t263 = cos(qJ(3));
t260 = sin(qJ(3));
t315 = qJD(1) * t260;
t223 = -t262 * t263 * qJD(1) + t259 * t315;
t219 = qJD(5) + t223;
t356 = t219 ^ 2;
t225 = (t259 * t263 + t260 * t262) * qJD(1);
t251 = qJD(3) + qJD(4);
t258 = sin(qJ(5));
t261 = cos(qJ(5));
t206 = t225 * t258 - t261 * t251;
t357 = t206 ^ 2;
t186 = t357 - t356;
t313 = qJD(1) * qJD(3);
t301 = t263 * t313;
t312 = t260 * qJDD(1);
t230 = t301 + t312;
t248 = t263 * qJDD(1);
t302 = t260 * t313;
t231 = t248 - t302;
t295 = t259 * t230 - t262 * t231;
t181 = -qJD(4) * t225 - t295;
t177 = qJDD(5) - t181;
t208 = t225 * t261 + t251 * t258;
t327 = t208 * t206;
t136 = -t327 - t177;
t338 = t136 * t258;
t105 = -t186 * t261 - t338;
t190 = t219 * t208;
t182 = -qJD(4) * t223 + t230 * t262 + t231 * t259;
t311 = qJDD(3) + qJDD(4);
t297 = -t258 * t182 + t261 * t311;
t284 = qJD(5) * t208 - t297;
t123 = -t190 + t284;
t420 = t260 * (t105 * t262 + t123 * t259) + t263 * (t105 * t259 - t123 * t262);
t205 = t208 ^ 2;
t368 = -t205 - t356;
t91 = t261 * t368 + t338;
t419 = pkin(2) * t91;
t418 = pkin(3) * t91;
t417 = pkin(4) * t91;
t416 = pkin(9) * t91;
t337 = t136 * t261;
t93 = -t258 * t368 + t337;
t415 = pkin(9) * t93;
t256 = cos(pkin(10));
t414 = t256 * t91;
t413 = t259 * t93;
t412 = t262 * t93;
t367 = t205 - t357;
t278 = -t261 * t182 - t258 * t311;
t274 = -t206 * qJD(5) - t278;
t328 = t206 * t219;
t363 = -t328 + t274;
t340 = t363 * t258;
t370 = t190 + t284;
t70 = t370 * t261 + t340;
t409 = t260 * (-t259 * t367 + t262 * t70) + t263 * (t259 * t70 + t262 * t367);
t255 = sin(pkin(10));
t364 = -t327 + t177;
t335 = t364 * t261;
t361 = -t356 - t357;
t372 = t258 * t361 + t335;
t336 = t364 * t258;
t371 = t261 * t361 - t336;
t388 = t259 * t370 + t262 * t371;
t389 = t259 * t371 - t262 * t370;
t402 = -t260 * t389 + t263 * t388;
t406 = pkin(1) * (t255 * t402 - t256 * t372) + pkin(7) * t402 - pkin(2) * t372;
t405 = pkin(3) * t389;
t404 = pkin(8) * t389;
t101 = -t186 * t258 + t337;
t403 = -pkin(3) * t372 + pkin(8) * t388;
t401 = t260 * t388 + t263 * t389;
t362 = t328 + t274;
t187 = -t205 + t356;
t390 = -t187 * t258 + t335;
t400 = t260 * (t259 * t362 + t262 * t390) + t263 * (t259 * t390 - t262 * t362);
t397 = pkin(4) * t372;
t396 = pkin(9) * t371;
t395 = pkin(9) * t372;
t391 = t261 * t187 + t336;
t339 = t363 * t261;
t68 = -t370 * t258 + t339;
t366 = t205 + t357;
t387 = pkin(4) * t366;
t386 = -qJ(6) * t258 - pkin(4);
t385 = qJ(6) * t363;
t382 = t259 * t366;
t201 = t225 * t223;
t369 = -t201 + t311;
t380 = t259 * t369;
t376 = t262 * t366;
t374 = t262 * t369;
t198 = pkin(4) * t223 - pkin(9) * t225;
t355 = t251 ^ 2;
t265 = qJD(1) ^ 2;
t353 = sin(qJ(1));
t354 = cos(qJ(1));
t279 = t353 * g(1) - t354 * g(2);
t275 = qJDD(1) * pkin(1) + t279;
t280 = t354 * g(1) + t353 * g(2);
t228 = -t265 * pkin(1) - t280;
t322 = t256 * t228;
t268 = -t265 * pkin(2) + qJDD(1) * pkin(7) + t255 * t275 + t322;
t317 = -g(3) + qJDD(2);
t179 = t260 * t317 + t263 * t268;
t237 = qJD(3) * pkin(3) - pkin(8) * t315;
t253 = t263 ^ 2;
t250 = t253 * t265;
t151 = -pkin(3) * t250 + t231 * pkin(8) - qJD(3) * t237 + t179;
t319 = t262 * t151;
t267 = t260 * t268;
t320 = t260 * t265;
t351 = t230 * pkin(8);
t365 = qJDD(3) * pkin(3) - t267 + (pkin(3) * t320 + pkin(8) * t313 + t317) * t263 - t351;
t97 = t259 * t365 + t319;
t79 = -t355 * pkin(4) + t311 * pkin(9) - t223 * t198 + t97;
t222 = t256 * t275;
t296 = -t255 * t228 + t222;
t193 = -qJDD(1) * pkin(2) - t265 * pkin(7) - t296;
t154 = -t231 * pkin(3) - pkin(8) * t250 + t237 * t315 + t193;
t216 = t251 * t223;
t161 = t182 - t216;
t82 = -t161 * pkin(9) + (t225 * t251 - t181) * pkin(4) + t154;
t45 = t258 * t79 - t261 * t82;
t46 = t258 * t82 + t261 * t79;
t20 = t258 * t45 + t261 * t46;
t167 = pkin(5) * t206 - qJ(6) * t208;
t294 = t177 * qJ(6) - t206 * t167 + t46;
t360 = -(t368 + t356) * pkin(5) - qJ(6) * t136 + t294;
t325 = t219 * t261;
t306 = t206 * t325;
t285 = t258 * t284 + t306;
t307 = t262 * t327;
t308 = t259 * t327;
t359 = t260 * (t262 * t285 - t308) + t263 * (t259 * t285 + t307);
t326 = t219 * t258;
t184 = t208 * t326;
t291 = t184 - t306;
t358 = t260 * (t177 * t259 + t262 * t291) + t263 * (-t262 * t177 + t259 * t291);
t220 = t223 ^ 2;
t221 = t225 ^ 2;
t352 = pkin(4) * t259;
t96 = t259 * t151 - t262 * t365;
t78 = -t311 * pkin(4) - t355 * pkin(9) + t225 * t198 + t96;
t350 = -pkin(4) * t78 + pkin(9) * t20;
t74 = t258 * t78;
t53 = t259 * t97 - t262 * t96;
t349 = t260 * t53;
t75 = t261 * t78;
t347 = qJ(6) * t261;
t342 = t362 * t258;
t341 = t362 * t261;
t334 = t154 * t259;
t333 = t154 * t262;
t196 = t201 + t311;
t330 = t196 * t259;
t329 = t196 * t262;
t324 = t251 * t259;
t323 = t251 * t262;
t241 = t263 * t320;
t235 = qJDD(3) + t241;
t321 = t260 * t235;
t236 = qJDD(3) - t241;
t318 = t263 * t236;
t314 = qJD(6) * t219;
t130 = (qJD(5) + t219) * t206 + t278;
t310 = pkin(4) * t130 + t415 + t74;
t309 = -pkin(4) * t370 + t396 - t75;
t305 = -pkin(1) * t256 - pkin(2);
t304 = pkin(1) * t255 + pkin(7);
t303 = -pkin(4) * t262 - pkin(3);
t54 = t259 * t96 + t262 * t97;
t214 = 0.2e1 * t314;
t290 = t214 + t294;
t24 = (t366 - t356) * pkin(5) + t290;
t39 = -t177 * pkin(5) - qJ(6) * t356 + t167 * t208 + qJDD(6) + t45;
t27 = qJ(6) * t366 + t39;
t71 = -t123 * t261 + t342;
t300 = pkin(9) * t71 + t261 * t24 + t258 * t27 + t387;
t125 = (-qJD(5) + t219) * t208 + t297;
t73 = t125 * t261 + t342;
t299 = pkin(9) * t73 + t20 + t387;
t178 = -t263 * t317 + t267;
t131 = t260 * t178 + t263 * t179;
t38 = -pkin(5) * t356 + t290;
t293 = -pkin(5) * t39 + qJ(6) * t38;
t292 = t206 * t326 - t261 * t284;
t289 = -pkin(5) * t362 - qJ(6) * t123;
t287 = t258 * t46 - t261 * t45;
t272 = t284 * pkin(5) - t385 + t78;
t271 = 0.2e1 * qJD(6) * t208 - t272;
t30 = -pkin(5) * t190 + t271 + t385;
t286 = pkin(4) * t363 + pkin(5) * t339 + t258 * t30 - t415;
t31 = (-t370 - t190) * pkin(5) + t271;
t283 = t261 * t31 + t386 * t370 + t396;
t281 = (-t206 * t258 - t208 * t261) * t219;
t277 = (-qJD(4) + t251) * t225 - t295;
t10 = t258 * t39 + t261 * t38;
t42 = (pkin(5) * t219 - 0.2e1 * qJD(6)) * t208 + t272;
t276 = pkin(9) * t10 + (-pkin(5) * t261 + t386) * t42;
t117 = t261 * t274 - t184;
t270 = t260 * (t262 * t117 + t308) + t263 * (t259 * t117 - t307);
t269 = pkin(5) * t364 + qJ(6) * t361 - t39;
t264 = qJD(3) ^ 2;
t252 = t260 ^ 2;
t249 = t252 * t265;
t239 = -t250 - t264;
t238 = -t249 - t264;
t234 = t249 + t250;
t233 = (t252 + t253) * qJDD(1);
t232 = t248 - 0.2e1 * t302;
t229 = 0.2e1 * t301 + t312;
t213 = -t221 + t355;
t212 = t220 - t355;
t211 = -t221 - t355;
t210 = -t238 * t260 - t318;
t209 = t239 * t263 - t321;
t200 = t221 - t220;
t194 = -t355 - t220;
t183 = -t220 - t221;
t164 = -t211 * t259 - t329;
t163 = t211 * t262 - t330;
t162 = t182 + t216;
t157 = (qJD(4) + t251) * t225 + t295;
t153 = t194 * t262 - t380;
t152 = t194 * t259 + t374;
t116 = t208 * t325 + t258 * t274;
t109 = -t163 * t260 + t164 * t263;
t108 = t162 * t259 + t262 * t277;
t107 = -t162 * t262 + t259 * t277;
t98 = -t152 * t260 + t153 * t263;
t69 = t125 * t258 - t341;
t67 = -t123 * t258 - t341;
t63 = -t107 * t260 + t108 * t263;
t61 = -t130 * t259 + t412;
t59 = t130 * t262 + t413;
t57 = -t259 * t363 - t412;
t55 = t262 * t363 - t413;
t52 = t262 * t73 - t382;
t51 = t262 * t71 - t382;
t50 = t259 * t73 + t376;
t49 = t259 * t71 + t376;
t48 = t75 - t416;
t47 = t74 - t395;
t40 = -pkin(4) * t67 - t289;
t37 = t46 - t417;
t35 = -t260 * t59 + t263 * t61;
t34 = t45 - t397;
t32 = -t260 * t55 + t263 * t57;
t25 = t263 * t54 - t349;
t21 = -t260 * t49 + t263 * t51;
t17 = -t269 - t397;
t16 = -t258 * t31 - t347 * t370 - t395;
t15 = -0.2e1 * t314 - t360 + t417;
t14 = -pkin(5) * t340 + t261 * t30 + t416;
t13 = t20 * t262 + t259 * t78;
t12 = t20 * t259 - t262 * t78;
t11 = -pkin(9) * t69 - t287;
t9 = t258 * t38 - t261 * t39;
t7 = -pkin(9) * t67 - t24 * t258 + t261 * t27;
t6 = t10 * t262 + t259 * t42;
t5 = t10 * t259 - t262 * t42;
t3 = -pkin(9) * t9 + (pkin(5) * t258 - t347) * t42;
t2 = -pkin(4) * t9 - t293;
t1 = -t260 * t5 + t263 * t6;
t4 = [0, 0, 0, 0, 0, qJDD(1), t279, t280, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t256 - t255 * t265) + t296, -t322 - t255 * t279 + (-0.2e1 * qJDD(1) * t255 - t256 * t265) * pkin(1), 0, pkin(1) * (t255 ^ 2 * t275 + t256 * t222), (t230 + t301) * t260, t229 * t263 + t232 * t260, t321 + t263 * (-t249 + t264), (t231 - t302) * t263, t260 * (t250 - t264) + t318, 0, -t263 * t193 + pkin(2) * t232 + pkin(7) * t209 + pkin(1) * (t209 * t255 + t232 * t256), t260 * t193 - pkin(2) * t229 + pkin(7) * t210 + pkin(1) * (t210 * t255 - t229 * t256), pkin(2) * t234 + pkin(7) * t233 + pkin(1) * (t233 * t255 + t234 * t256) + t131, -pkin(2) * t193 + pkin(7) * t131 + pkin(1) * (t131 * t255 - t193 * t256), t260 * (t182 * t262 - t225 * t324) + t263 * (t182 * t259 + t225 * t323), t260 * (-t157 * t262 - t161 * t259) + t263 * (-t157 * t259 + t161 * t262), t260 * (-t213 * t259 + t374) + t263 * (t213 * t262 + t380), t260 * (-t181 * t259 + t223 * t323) + t263 * (t181 * t262 + t223 * t324), t260 * (t212 * t262 - t330) + t263 * (t212 * t259 + t329), (t260 * (-t223 * t262 + t225 * t259) + t263 * (-t223 * t259 - t225 * t262)) * t251, t260 * (-pkin(8) * t152 + t334) + t263 * (-pkin(3) * t157 + pkin(8) * t153 - t333) - pkin(2) * t157 + pkin(7) * t98 + pkin(1) * (-t157 * t256 + t255 * t98), t260 * (-pkin(8) * t163 + t333) + t263 * (-pkin(3) * t161 + pkin(8) * t164 + t334) - pkin(2) * t161 + pkin(7) * t109 + pkin(1) * (t109 * t255 - t161 * t256), t260 * (-pkin(8) * t107 - t53) + t263 * (-pkin(3) * t183 + pkin(8) * t108 + t54) - pkin(2) * t183 + pkin(7) * t63 + pkin(1) * (-t183 * t256 + t255 * t63), -pkin(8) * t349 + t263 * (-pkin(3) * t154 + pkin(8) * t54) - pkin(2) * t154 + pkin(7) * t25 + pkin(1) * (-t154 * t256 + t25 * t255), t270, -t409, t400, t359, -t420, t358, t260 * (-t259 * t34 + t262 * t47 - t404) + t263 * (t259 * t47 + t262 * t34 + t403) + t406, t260 * (-pkin(8) * t59 - t259 * t37 + t262 * t48) + t263 * (pkin(8) * t61 + t259 * t48 + t262 * t37 - t418) - t419 + pkin(7) * t35 + pkin(1) * (t255 * t35 - t414), t260 * (-pkin(8) * t50 + t11 * t262) + t263 * (pkin(8) * t52 + t259 * t11) + t304 * (-t260 * t50 + t263 * t52) + (t260 * t352 + t263 * t303 + t305) * t69, (t260 * (-pkin(9) * t262 + t352) + t263 * (-pkin(9) * t259 + t303) + t305) * t287 + (t304 + pkin(8)) * (-t12 * t260 + t13 * t263), t270, t400, t409, t358, t420, t359, t260 * (t16 * t262 - t17 * t259 - t404) + t263 * (t16 * t259 + t17 * t262 + t403) + t406, t260 * (-pkin(8) * t49 - t259 * t40 + t262 * t7) + t263 * (-pkin(3) * t67 + pkin(8) * t51 + t259 * t7 + t262 * t40) - pkin(2) * t67 + pkin(7) * t21 + pkin(1) * (t21 * t255 - t256 * t67), t260 * (-pkin(8) * t55 + t14 * t262 - t15 * t259) + t263 * (pkin(8) * t57 + t14 * t259 + t15 * t262 + t418) + t419 + pkin(7) * t32 + pkin(1) * (t255 * t32 + t414), t260 * (-pkin(8) * t5 - t2 * t259 + t262 * t3) + t263 * (-pkin(3) * t9 + pkin(8) * t6 + t2 * t262 + t259 * t3) - pkin(2) * t9 + pkin(7) * t1 + pkin(1) * (t1 * t255 - t256 * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, 0, 0, 0, 0, 0, 0, t235 * t263 + t239 * t260, -t236 * t260 + t238 * t263, 0, -t178 * t263 + t179 * t260, 0, 0, 0, 0, 0, 0, t152 * t263 + t153 * t260, t163 * t263 + t164 * t260, t107 * t263 + t108 * t260, t260 * t54 + t263 * t53, 0, 0, 0, 0, 0, 0, t401, t260 * t61 + t263 * t59, t260 * t52 + t263 * t50, t12 * t263 + t13 * t260, 0, 0, 0, 0, 0, 0, t401, t260 * t51 + t263 * t49, t260 * t57 + t263 * t55, t260 * t6 + t263 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, t249 - t250, t312, t241, t248, qJDD(3), -t178, -t179, 0, 0, t201, t200, t162, -t201, t277, t311, pkin(3) * t152 - t96, -t319 - t259 * (pkin(8) * t301 - t178 - t351) + (-t235 * t259 + t163) * pkin(3), pkin(3) * t107, pkin(3) * t53, t116, t68, t391, t292, -t101, t281, t309 + t405, pkin(3) * t59 + t310, pkin(3) * t50 + t299, pkin(3) * t12 + t350, t116, t391, -t68, t281, t101, t292, t283 + t405, pkin(3) * t49 + t300, pkin(3) * t55 + t286, pkin(3) * t5 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t200, t162, -t201, t277, t311, -t96, -t97, 0, 0, t116, t68, t391, t292, -t101, t281, t309, t310, t299, t350, t116, t391, -t68, t281, t101, t292, t283, t300, t286, t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t327, t367, t362, -t327, -t123, t177, -t45, -t46, 0, 0, t327, t362, -t367, t177, t123, -t327, t269, t289, t214 + t360, t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t364, t362, t368, t39;];
tauJ_reg  = t4;
