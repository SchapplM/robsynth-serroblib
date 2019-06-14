% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 08:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRPRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:35:54
% EndTime: 2019-05-07 08:36:15
% DurationCPUTime: 6.93s
% Computational Cost: add. (23588->459), mult. (47051->540), div. (0->0), fcn. (32215->8), ass. (0->304)
t268 = cos(qJ(3));
t264 = sin(qJ(3));
t265 = sin(qJ(2));
t326 = qJD(1) * qJD(2);
t251 = t265 * t326;
t269 = cos(qJ(2));
t325 = t269 * qJDD(1);
t238 = -t251 + t325;
t231 = -qJDD(3) + t238;
t331 = qJD(1) * t265;
t232 = -t268 * qJD(2) + t264 * t331;
t234 = t264 * qJD(2) + t268 * t331;
t340 = t234 * t232;
t280 = t231 - t340;
t346 = t280 * t264;
t230 = t234 ^ 2;
t249 = qJD(1) * t269 - qJD(3);
t383 = t249 ^ 2;
t396 = -t230 - t383;
t140 = t268 * t396 + t346;
t430 = pkin(1) * t140;
t429 = pkin(2) * t140;
t428 = pkin(8) * t140;
t345 = t280 * t268;
t142 = -t264 * t396 + t345;
t427 = pkin(8) * t142;
t426 = t142 * t269;
t263 = sin(qJ(5));
t267 = cos(qJ(5));
t196 = -t267 * t232 + t234 * t263;
t198 = t232 * t263 + t234 * t267;
t146 = t198 * t196;
t225 = qJDD(5) + t231;
t405 = t146 - t225;
t130 = pkin(5) * t405;
t217 = t234 * t249;
t252 = t265 * qJDD(1);
t320 = t269 * t326;
t237 = t252 + t320;
t313 = -t268 * qJDD(2) + t264 * t237;
t285 = qJD(3) * t234 + t313;
t154 = t217 + t285;
t384 = t232 ^ 2;
t210 = t384 - t383;
t425 = t269 * t154 + t265 * (t210 * t268 + t346);
t291 = -t264 * qJDD(2) - t268 * t237;
t190 = -qJD(3) * t232 - t291;
t341 = t232 * t249;
t398 = t190 + t341;
t351 = t398 * t264;
t395 = t230 - t384;
t397 = -t217 + t285;
t424 = t265 * (t397 * t268 + t351) + t269 * t395;
t381 = pkin(3) + pkin(4);
t177 = t231 + t340;
t343 = t177 * t268;
t392 = -t383 - t384;
t404 = t264 * t392 - t343;
t423 = pkin(2) * t404;
t344 = t177 * t264;
t403 = t268 * t392 + t344;
t422 = pkin(8) * t403;
t421 = pkin(8) * t404;
t420 = qJ(4) * t398;
t211 = -t230 + t383;
t419 = t268 * t211 - t344;
t417 = t210 * t264 - t345;
t399 = t190 - t341;
t415 = t265 * (-t211 * t264 - t343) - t269 * t399;
t414 = pkin(7) * (t265 * t397 + t269 * t403) - pkin(1) * t404;
t394 = t230 + t384;
t413 = pkin(2) * t394;
t122 = -t196 * qJD(5) + t267 * t190 + t263 * t285;
t245 = qJD(5) + t249;
t172 = t245 * t196;
t104 = t172 + t122;
t412 = qJ(6) * t104;
t410 = t265 * t394;
t355 = t405 * t263;
t354 = t405 * t267;
t271 = qJD(1) ^ 2;
t266 = sin(qJ(1));
t270 = cos(qJ(1));
t306 = g(1) * t270 + g(2) * t266;
t358 = qJDD(1) * pkin(7);
t223 = -pkin(1) * t271 - t306 + t358;
t309 = -pkin(2) * t269 - pkin(8) * t265;
t312 = t271 * t309 + t223;
t371 = t269 * g(3);
t382 = qJD(2) ^ 2;
t165 = -qJDD(2) * pkin(2) - t382 * pkin(8) + t312 * t265 + t371;
t277 = t285 * pkin(3) + t165 - t420;
t402 = -t397 * t264 + t268 * t398;
t400 = -t172 + t122;
t307 = pkin(4) * t249 - pkin(9) * t234;
t393 = t384 * pkin(9) - t234 * t307;
t167 = pkin(5) * t245 - qJ(6) * t198;
t391 = -t198 * t167 - qJDD(6);
t194 = t196 ^ 2;
t243 = t245 ^ 2;
t139 = -t243 - t194;
t91 = t139 * t263 - t354;
t92 = t139 * t267 + t355;
t390 = qJ(4) * t92 - t381 * t91;
t195 = t198 ^ 2;
t163 = -t195 - t243;
t132 = t146 + t225;
t357 = t132 * t263;
t109 = t163 * t267 - t357;
t356 = t132 * t267;
t110 = -t163 * t263 - t356;
t389 = qJ(4) * t110 - t109 * t381;
t318 = g(1) * t266 - t270 * g(2);
t222 = qJDD(1) * pkin(1) + pkin(7) * t271 + t318;
t295 = -t238 + t251;
t296 = t237 + t320;
t151 = t295 * pkin(2) - t296 * pkin(8) - t222;
t374 = g(3) * t265;
t166 = -t382 * pkin(2) + qJDD(2) * pkin(8) + t312 * t269 - t374;
t125 = t264 * t151 + t268 * t166;
t199 = pkin(3) * t232 - qJ(4) * t234;
t311 = t231 * qJ(4) + t232 * t199 - t125;
t387 = -pkin(3) * (t396 + t383) - qJ(4) * t280 - t311;
t329 = qJD(6) * t196;
t186 = -0.2e1 * t329;
t315 = t263 * t190 - t267 * t285;
t121 = -qJD(5) * t198 - t315;
t124 = -t268 * t151 + t166 * t264;
t87 = t231 * pkin(3) - qJ(4) * t383 + t199 * t234 + qJDD(4) + t124;
t70 = t177 * pkin(4) - pkin(9) * t399 + t87;
t330 = qJD(4) * t249;
t242 = -0.2e1 * t330;
t294 = t242 - t311;
t85 = -pkin(3) * t383 + t294;
t74 = -pkin(4) * t384 + t285 * pkin(9) - t249 * t307 + t85;
t40 = t263 * t70 + t267 * t74;
t302 = t194 * pkin(5) - t121 * qJ(6) + t245 * t167 - t40;
t26 = t186 - t302;
t328 = qJD(6) * t198;
t188 = -0.2e1 * t328;
t39 = t263 * t74 - t267 * t70;
t282 = t39 + t412 + t130;
t25 = t188 - t282;
t363 = t25 * t267;
t11 = t26 * t263 + t363;
t364 = t25 * t263;
t12 = t26 * t267 - t364;
t24 = pkin(5) * t25;
t386 = qJ(4) * t12 - t381 * t11 - t24;
t19 = t263 * t40 - t267 * t39;
t20 = t263 * t39 + t267 * t40;
t385 = qJ(4) * t20 - t381 * t19;
t101 = (qJD(5) - t245) * t198 + t315;
t53 = t264 * t92 - t268 * t91;
t380 = pkin(2) * t53;
t66 = -t109 * t268 + t110 * t264;
t379 = pkin(2) * t66;
t61 = -t101 * t263 - t104 * t267;
t63 = -t101 * t267 + t104 * t263;
t35 = t264 * t63 - t268 * t61;
t378 = pkin(8) * t35;
t377 = pkin(8) * t53;
t376 = pkin(8) * t66;
t375 = pkin(3) * t268;
t373 = t121 * pkin(5);
t127 = -t194 - t195;
t36 = t264 * t61 + t268 * t63;
t370 = pkin(7) * (-t127 * t265 + t269 * t36) - pkin(1) * t35;
t100 = (qJD(5) + t245) * t198 + t315;
t54 = t264 * t91 + t268 * t92;
t369 = pkin(7) * (-t100 * t265 + t269 * t54) - pkin(1) * t53;
t67 = t109 * t264 + t110 * t268;
t368 = pkin(7) * (-t265 * t400 + t269 * t67) - pkin(1) * t66;
t367 = pkin(2) * t100 + pkin(8) * t54;
t366 = pkin(2) * t400 + pkin(8) * t67;
t316 = -pkin(3) * t249 - 0.2e1 * qJD(4);
t75 = t313 * pkin(4) + (pkin(4) * qJD(3) + t316) * t234 + t277 + t393;
t362 = t263 * t75;
t361 = t267 * t75;
t360 = pkin(2) * t127 + pkin(8) * t36;
t359 = qJ(4) * t268;
t350 = t399 * t264;
t349 = t399 * t268;
t348 = t165 * t264;
t347 = t165 * t268;
t339 = t245 * t263;
t338 = t245 * t267;
t337 = t249 * t264;
t336 = t249 * t268;
t248 = t269 * t271 * t265;
t335 = t265 * (qJDD(2) + t248);
t333 = t269 * (-t248 + qJDD(2));
t324 = qJ(4) * t63 - t381 * t61;
t323 = t232 * t336;
t322 = t269 * t146;
t321 = t269 * t340;
t319 = -qJ(4) * t264 - pkin(2);
t81 = t124 * t264 + t268 * t125;
t207 = t265 * t223 + t371;
t208 = t223 * t269 - t374;
t314 = t265 * t207 + t269 * t208;
t308 = -pkin(3) * t87 + qJ(4) * t85;
t209 = t234 * t337;
t305 = t265 * (t190 * t268 + t209) - t321;
t304 = -t232 * t337 - t268 * t285;
t303 = t40 + t389;
t301 = -pkin(9) * t91 + qJ(4) * t100;
t300 = -pkin(9) * t61 + qJ(4) * t127;
t299 = -pkin(2) * t35 - t324;
t298 = -pkin(3) * t399 - qJ(4) * t154;
t297 = -pkin(9) * t109 + qJ(4) * t400;
t293 = t124 * t268 - t125 * t264;
t292 = -pkin(1) + t309;
t290 = -pkin(9) * t92 + t381 * t100;
t289 = -pkin(9) * t63 + t381 * t127;
t288 = pkin(5) * t163 + t302;
t287 = -pkin(9) * t110 + t381 * t400;
t286 = t39 + t390;
t283 = (t232 * t264 + t234 * t268) * t249;
t281 = t265 * (-t209 + t323) + t269 * t231;
t279 = -t288 + t389;
t278 = t130 + t282 + t390;
t276 = 0.2e1 * qJD(4) * t234 - t277;
t275 = -pkin(3) * t177 + qJ(4) * t392 - t87;
t274 = pkin(3) * t217 + t276;
t273 = t265 * (t264 * t285 - t323) + t321;
t272 = t75 + t391;
t260 = t269 ^ 2;
t259 = t265 ^ 2;
t257 = t260 * t271;
t255 = t259 * t271;
t239 = -0.2e1 * t251 + t325;
t236 = t252 + 0.2e1 * t320;
t189 = 0.2e1 * t328;
t187 = 0.2e1 * t329;
t169 = -t195 + t243;
t168 = t194 - t243;
t160 = (qJD(3) - t249) * t232 + t291;
t155 = (-qJD(3) - t249) * t234 - t313;
t148 = t190 * t264 - t234 * t336;
t144 = t195 - t194;
t129 = (-t196 * t267 + t198 * t263) * t245;
t128 = (t196 * t263 + t198 * t267) * t245;
t118 = t155 * t268 + t350;
t117 = -t154 * t268 + t350;
t115 = -t154 * t264 - t349;
t114 = t168 * t267 - t357;
t113 = -t169 * t263 - t354;
t112 = -t168 * t263 - t356;
t111 = -t169 * t267 + t355;
t97 = pkin(5) * t104;
t96 = t122 * t267 - t198 * t339;
t95 = -t122 * t263 - t198 * t338;
t94 = -t121 * t263 + t196 * t338;
t93 = -t121 * t267 - t196 * t339;
t86 = t316 * t234 + t277;
t84 = t128 * t268 + t129 * t264;
t83 = t265 * (-t128 * t264 + t129 * t268) + t269 * t225;
t82 = qJ(4) * t394 + t87;
t79 = (t394 - t383) * pkin(3) + t294;
t78 = (-t397 + t217) * pkin(3) + t276;
t77 = t274 + t420;
t76 = -pkin(5) * t400 - qJ(6) * t132;
t72 = t112 * t268 + t114 * t264;
t71 = t111 * t268 + t113 * t264;
t62 = -t100 * t267 - t263 * t400;
t60 = t100 * t263 - t267 * t400;
t56 = t264 * t96 + t268 * t95;
t55 = t264 * t94 + t268 * t93;
t50 = t264 * t87 + t268 * t85;
t49 = t264 * t85 - t268 * t87;
t48 = t265 * (-t264 * t95 + t268 * t96) + t322;
t47 = t265 * (-t264 * t93 + t268 * t94) - t322;
t46 = t265 * (-t112 * t264 + t114 * t268) - t269 * t101;
t45 = t265 * (-t111 * t264 + t113 * t268) + t269 * t104;
t43 = t194 * qJ(6) + t272 + t373;
t41 = -t285 * pkin(4) + (-t194 - t163) * qJ(6) + t274 - t373 - t391 - t393;
t38 = t297 - t361;
t37 = t301 - t362;
t34 = t264 * t62 + t268 * t60;
t31 = t265 * (-t264 * t60 + t268 * t62) + t269 * t144;
t30 = t287 + t362;
t29 = (t194 + t139) * qJ(6) + (t121 - t100) * pkin(5) + t272;
t28 = t290 - t361;
t23 = t189 + t282 + t412;
t22 = -pkin(5) * t127 - qJ(6) * t101 + t26;
t21 = -t263 * t76 + t267 * t41 + t297;
t18 = qJ(6) * t354 - t263 * t29 + t301;
t17 = -t263 * t41 - t267 * t76 + t287;
t16 = pkin(5) * t43 + qJ(6) * t26;
t15 = -qJ(6) * t355 - t267 * t29 + t290;
t14 = -pkin(9) * t19 - qJ(4) * t75;
t13 = -t19 + t300;
t10 = -t20 + t289;
t9 = -pkin(9) * t20 - t381 * t75;
t8 = t19 * t264 + t20 * t268;
t7 = -t19 * t268 + t20 * t264;
t6 = -t22 * t263 + t23 * t267 + t300;
t5 = -t22 * t267 - t23 * t263 + t289;
t4 = t11 * t264 + t12 * t268;
t3 = -t11 * t268 + t12 * t264;
t2 = -pkin(9) * t11 - qJ(4) * t43 - qJ(6) * t363 - t16 * t263;
t1 = -pkin(9) * t12 + qJ(6) * t364 - t16 * t267 - t381 * t43;
t27 = [0, 0, 0, 0, 0, qJDD(1), t318, t306, 0, 0, t296 * t265, t236 * t269 + t239 * t265, t335 + t269 * (-t255 + t382), -t295 * t269, t265 * (t257 - t382) + t333, 0, t269 * t222 + pkin(1) * t239 + pkin(7) * (t269 * (-t257 - t382) - t335), -t265 * t222 - pkin(1) * t236 + pkin(7) * (-t333 - t265 * (-t255 - t382)), pkin(1) * (t255 + t257) + (t259 + t260) * t358 + t314, pkin(1) * t222 + pkin(7) * t314, t305, -t424, t415, t273, t425, t281, t265 * (t348 - t421) + t269 * (t124 - t423) + t414, t265 * (t347 - t428) + t269 * (t125 - t429) - t430 + pkin(7) * (-t160 * t265 + t426), t265 * t293 + pkin(7) * (t118 * t269 - t410) + t292 * (t155 * t264 - t349), pkin(7) * (t165 * t265 + t269 * t81) - t292 * t293, t305, t415, t424, t281, -t425, t273, t265 * (-t264 * t78 - t359 * t397 - t421) + t269 * (-t275 - t423) + t414, t265 * (-pkin(8) * t115 - t264 * t79 + t268 * t82) + t269 * (-pkin(2) * t115 - t298) - pkin(1) * t115 + pkin(7) * (t117 * t269 - t410), t265 * (-pkin(3) * t351 + t268 * t77 + t428) + t269 * (0.2e1 * t330 - t387 + t429) + t430 + pkin(7) * (-t265 * t398 - t426), t265 * (-pkin(8) * t49 + (pkin(3) * t264 - t359) * t86) + t269 * (-pkin(2) * t49 - t308) - pkin(1) * t49 + pkin(7) * (t265 * t86 + t269 * t50), t48, t31, t45, t47, t46, t83, t265 * (-t264 * t28 + t268 * t37 - t377) + t269 * (-t286 - t380) + t369, t265 * (-t264 * t30 + t268 * t38 - t376) + t269 * (-t303 - t379) + t368, t265 * (-t10 * t264 + t13 * t268 - t378) + t269 * t299 + t370, t265 * (-pkin(8) * t7 + t14 * t268 - t264 * t9) + t269 * (-pkin(2) * t7 - t385) - pkin(1) * t7 + pkin(7) * (t265 * t75 + t269 * t8), t48, t31, t45, t47, t46, t83, t265 * (-t15 * t264 + t18 * t268 - t377) + t269 * (t188 - t278 - t380) + t369, t265 * (-t17 * t264 + t21 * t268 - t376) + t269 * (t187 - t279 - t379) + t368, t265 * (-t264 * t5 + t268 * t6 - t378) + t269 * (t299 - t97) + t370, t265 * (-pkin(8) * t3 - t1 * t264 + t2 * t268) + t269 * (-pkin(2) * t3 - t386) - pkin(1) * t3 + pkin(7) * (t265 * t43 + t269 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t248, t255 - t257, t252, t248, t325, qJDD(2), -t207, -t208, 0, 0, t148, t402, t419, t304, t417, t283, -pkin(2) * t397 - t347 + t422, pkin(2) * t160 + t348 + t427, pkin(8) * t118 + t413 + t81, -pkin(2) * t165 + pkin(8) * t81, t148, t419, -t402, t283, -t417, t304, t268 * t78 + t319 * t397 + t422, pkin(8) * t117 + t264 * t82 + t268 * t79 + t413, -t427 + t264 * t77 + (pkin(2) + t375) * t398, pkin(8) * t50 + (t319 - t375) * t86, t56, t34, t71, t55, t72, t84, t264 * t37 + t268 * t28 + t367, t264 * t38 + t268 * t30 + t366, t10 * t268 + t13 * t264 + t360, -pkin(2) * t75 + pkin(8) * t8 + t14 * t264 + t268 * t9, t56, t34, t71, t55, t72, t84, t15 * t268 + t18 * t264 + t367, t17 * t268 + t21 * t264 + t366, t264 * t6 + t268 * t5 + t360, -pkin(2) * t43 + pkin(8) * t4 + t1 * t268 + t2 * t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, t395, t399, -t340, -t154, -t231, -t124, -t125, 0, 0, t340, t399, -t395, -t231, t154, -t340, t275, t298, t242 + t387, t308, -t146, -t144, -t104, t146, t101, -t225, t286, t303, t324, t385, -t146, -t144, -t104, t146, t101, -t225, t189 + t278, t186 + t279, t97 + t324, t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t399, t396, t87, 0, 0, 0, 0, 0, 0, t91, t109, t61, t19, 0, 0, 0, 0, 0, 0, t91, t109, t61, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t144, t104, -t146, -t101, t225, -t39, -t40, 0, 0, t146, t144, t104, -t146, -t101, t225, -t130 + t25, t187 + t288, -t97, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t400, t127, -t43;];
tauJ_reg  = t27;
