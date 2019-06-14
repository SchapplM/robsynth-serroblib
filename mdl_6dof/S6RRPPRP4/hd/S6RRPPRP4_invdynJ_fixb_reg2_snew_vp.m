% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPPRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:21:52
% EndTime: 2019-05-06 09:22:18
% DurationCPUTime: 8.99s
% Computational Cost: add. (17425->443), mult. (38624->527), div. (0->0), fcn. (25968->8), ass. (0->288)
t261 = sin(qJ(2));
t264 = cos(qJ(2));
t258 = sin(pkin(9));
t259 = cos(pkin(9));
t323 = qJD(1) * t261;
t227 = -t259 * qJD(2) + t258 * t323;
t229 = qJD(2) * t258 + t259 * t323;
t260 = sin(qJ(5));
t263 = cos(qJ(5));
t188 = -t263 * t227 + t229 * t260;
t249 = t261 * qJDD(1);
t317 = qJD(1) * qJD(2);
t309 = t264 * t317;
t234 = t249 + t309;
t203 = qJDD(2) * t258 + t234 * t259;
t302 = -t259 * qJDD(2) + t234 * t258;
t122 = -t188 * qJD(5) + t263 * t203 + t260 * t302;
t320 = t264 * qJD(1);
t243 = qJD(5) + t320;
t337 = t188 * t243;
t391 = t122 - t337;
t245 = t261 * t317;
t316 = t264 * qJDD(1);
t235 = -t245 + t316;
t230 = qJDD(5) + t235;
t190 = t227 * t260 + t229 * t263;
t336 = t190 * t188;
t129 = -t336 - t230;
t355 = t129 * t260;
t187 = t190 ^ 2;
t374 = t243 ^ 2;
t390 = -t187 - t374;
t78 = t263 * t390 + t355;
t354 = t129 * t263;
t80 = -t260 * t390 + t354;
t46 = t258 * t80 - t259 * t78;
t48 = t258 * t78 + t259 * t80;
t478 = -pkin(7) * (-t261 * t391 + t264 * t48) + pkin(1) * t46;
t476 = qJ(3) * t46;
t370 = pkin(3) + pkin(4);
t475 = pkin(2) * t46 + qJ(4) * t80 - t370 * t78;
t474 = -pkin(2) * t391 - qJ(3) * t48;
t303 = t260 * t203 - t263 * t302;
t101 = (qJD(5) - t243) * t190 + t303;
t376 = t188 ^ 2;
t161 = t376 - t374;
t91 = t161 * t260 - t354;
t95 = -t161 * t263 - t355;
t472 = -t264 * t101 + t261 * (t258 * t91 - t259 * t95);
t471 = pkin(8) * t78;
t470 = pkin(8) * t80;
t334 = t227 * t229;
t279 = t235 - t334;
t341 = t279 * t258;
t225 = t229 ^ 2;
t256 = t264 ^ 2;
t266 = qJD(1) ^ 2;
t251 = t256 * t266;
t389 = -t225 - t251;
t143 = t259 * t389 + t341;
t340 = t279 * t259;
t145 = -t258 * t389 + t340;
t310 = t227 * t320;
t274 = t203 + t310;
t469 = pkin(7) * (t145 * t264 + t261 * t274) - pkin(1) * t143;
t467 = t258 * t95 + t259 * t91;
t464 = pkin(2) * t143;
t463 = qJ(3) * t143;
t462 = qJ(3) * t145;
t385 = t337 + t122;
t419 = -t101 * t263 + t385 * t260;
t420 = -t101 * t260 - t385 * t263;
t440 = t258 * t419 - t259 * t420;
t461 = qJ(3) * t440;
t460 = -pkin(2) * t440 - qJ(4) * t419 + t370 * t420;
t111 = -t376 - t187;
t441 = t258 * t420 + t259 * t419;
t459 = pkin(2) * t111 + qJ(3) * t441;
t138 = t187 - t376;
t100 = (qJD(5) + t243) * t190 + t303;
t55 = -t100 * t260 + t263 * t391;
t59 = t100 * t263 + t260 * t391;
t457 = t264 * t138 + t261 * (t258 * t55 - t259 * t59);
t456 = pkin(7) * (-t111 * t261 + t264 * t441) - pkin(1) * t440;
t386 = -t336 + t230;
t353 = t386 * t260;
t384 = -t374 - t376;
t394 = t263 * t384 - t353;
t352 = t386 * t263;
t395 = t260 * t384 + t352;
t411 = t258 * t394 - t259 * t395;
t455 = pkin(1) * t411;
t211 = t229 * t320;
t399 = t211 - t302;
t454 = pkin(3) * t399;
t453 = qJ(3) * t411;
t412 = t258 * t395 + t259 * t394;
t452 = qJ(3) * t412;
t451 = t264 * t412;
t450 = -pkin(2) * t411 - qJ(4) * t394 + t370 * t395;
t449 = -pkin(8) * t419 + t370 * t111;
t448 = -pkin(8) * t420 + qJ(4) * t111;
t375 = t227 ^ 2;
t165 = -t375 - t225;
t174 = -t302 - t211;
t275 = -t203 + t310;
t418 = t174 * t259 - t275 * t258;
t447 = -pkin(2) * t165 + qJ(3) * t418;
t207 = t375 - t251;
t446 = -t264 * t174 + t261 * (t207 * t259 + t341);
t445 = pkin(7) * (t165 * t261 + t264 * t418);
t444 = t258 * t59 + t259 * t55;
t162 = -t187 + t374;
t422 = -t162 * t260 + t352;
t423 = -t162 * t263 - t353;
t439 = t258 * t422 + t259 * t423;
t438 = t261 * (-t258 * t423 + t259 * t422) + t264 * t385;
t179 = t235 + t334;
t338 = t179 * t259;
t383 = -t375 - t251;
t397 = t258 * t383 - t338;
t435 = pkin(2) * t397;
t434 = pkin(8) * t394;
t433 = pkin(8) * t395;
t432 = qJ(3) * t397;
t208 = -t225 + t251;
t339 = t179 * t258;
t424 = t259 * t208 - t339;
t396 = t259 * t383 + t339;
t421 = pkin(2) * t399 + qJ(3) * t396;
t117 = t174 * t258 + t275 * t259;
t417 = t207 * t258 - t340;
t342 = t274 * t258;
t346 = t399 * t259;
t415 = t261 * (t342 - t346) + t264 * (t225 - t375);
t414 = t261 * (-t208 * t258 - t338) + t264 * t275;
t410 = pkin(7) * (-t261 * t399 + t264 * t396) - pkin(1) * t397;
t400 = t391 * qJ(6);
t262 = sin(qJ(1));
t265 = cos(qJ(1));
t296 = g(1) * t265 + g(2) * t262;
t360 = qJDD(1) * pkin(7);
t218 = -pkin(1) * t266 - t296 + t360;
t291 = -pkin(2) * t264 - qJ(3) * t261;
t232 = t291 * qJD(1);
t301 = qJD(1) * t232 + t218;
t392 = t301 * t261;
t193 = pkin(3) * t227 - qJ(4) * t229;
t307 = t262 * g(1) - t265 * g(2);
t217 = qJDD(1) * pkin(1) + t266 * pkin(7) + t307;
t290 = t234 + t309;
t149 = -t290 * qJ(3) + (-t235 + t245) * pkin(2) - t217;
t366 = t261 * g(3);
t373 = qJD(2) ^ 2;
t156 = -t373 * pkin(2) + qJDD(2) * qJ(3) + t264 * t301 - t366;
t324 = t258 * t149 + t259 * t156;
t387 = -t235 * qJ(4) - 0.2e1 * qJD(4) * t320 - t227 * t193 + t324;
t347 = t399 * t258;
t382 = t259 * t274 + t347;
t365 = t264 * g(3);
t287 = -qJDD(2) * pkin(2) - t373 * qJ(3) + qJDD(3) + t365;
t277 = t203 * qJ(4) - t287 + t454;
t329 = t261 * t218;
t333 = t227 * t264;
t381 = -qJD(1) * (qJ(4) * t333 - t232 * t261) - t277 + t329;
t121 = -qJD(5) * t190 - t303;
t331 = t243 * t260;
t315 = t188 * t331;
t280 = -t121 * t263 - t315;
t330 = t243 * t263;
t314 = t188 * t330;
t281 = -t121 * t260 + t314;
t380 = t258 * t281 + t259 * t280;
t160 = t190 * t330;
t293 = t160 + t315;
t159 = t190 * t331;
t294 = t159 - t314;
t379 = t258 * t294 + t259 * t293;
t313 = t264 * t336;
t378 = t261 * (-t258 * t280 + t259 * t281) - t313;
t377 = t261 * (-t258 * t293 + t259 * t294) + t264 * t230;
t372 = 2 * qJD(3);
t371 = 2 * qJD(6);
t369 = pkin(5) * t260;
t368 = pkin(5) * t263;
t367 = t121 * pkin(5);
t137 = pkin(5) * t188 - qJ(6) * t190;
t305 = -t259 * t149 + t258 * t156;
t278 = t235 * pkin(3) - qJ(4) * t251 + qJDD(4) + t305;
t319 = t372 + t193;
t269 = t235 * pkin(4) + t275 * pkin(8) + (pkin(4) * t227 + t319) * t229 + t278;
t286 = pkin(4) * t320 - pkin(8) * t229;
t322 = qJD(3) * t227;
t215 = -0.2e1 * t322;
t283 = t215 + t387;
t69 = -pkin(3) * t251 + t283;
t65 = -pkin(4) * t375 + pkin(8) * t302 - t286 * t320 + t69;
t32 = t260 * t269 + t263 * t65;
t292 = t230 * qJ(6) - t188 * t137 + t243 * t371 + t32;
t24 = -pkin(5) * t374 + t292;
t31 = t260 * t65 - t263 * t269;
t26 = -t230 * pkin(5) - qJ(6) * t374 + t190 * t137 + qJDD(6) + t31;
t364 = -pkin(5) * t26 + qJ(6) * t24;
t363 = -pkin(5) * t385 - qJ(6) * t101;
t321 = qJD(4) * t229;
t77 = -0.2e1 * t321 + t381;
t68 = pkin(4) * t302 + pkin(8) * t375 - t229 * t286 + t77;
t362 = t260 * t68;
t361 = t263 * t68;
t155 = t287 + t392;
t351 = t155 * t258;
t350 = t155 * t259;
t332 = t243 * t190;
t242 = t264 * t266 * t261;
t328 = t261 * (qJDD(2) + t242);
t325 = t264 * (-t242 + qJDD(2));
t312 = t229 * t333;
t109 = t215 + t324;
t311 = pkin(3) * t259 + pkin(2);
t108 = t229 * t372 + t305;
t64 = t108 * t258 + t259 * t109;
t306 = -qJ(6) * t263 + qJ(4);
t198 = t329 + t365;
t199 = t218 * t264 - t366;
t304 = t261 * t198 + t264 * t199;
t299 = qJ(6) * t260 + t370;
t298 = t258 * t310;
t86 = -t122 * t260 - t160;
t87 = t122 * t263 - t159;
t297 = t261 * (-t258 * t86 + t259 * t87) + t313;
t295 = t261 * (t203 * t259 + t258 * t211) - t312;
t17 = t260 * t32 - t263 * t31;
t18 = t260 * t31 + t263 * t32;
t289 = t108 * t259 - t109 * t258;
t288 = -pkin(1) + t291;
t285 = -t259 * t302 - t298;
t204 = t259 * t211;
t284 = t204 + t298;
t276 = -pkin(5) * t390 - qJ(6) * t129 + t24;
t220 = t264 * t235;
t273 = t220 + t261 * (t227 * t259 - t229 * t258) * t320;
t70 = t319 * t229 + t278;
t271 = pkin(5) * t386 + qJ(6) * t384 - t26;
t270 = t261 * (t258 * t302 - t259 * t310) + t312;
t268 = -pkin(5) * t332 + t190 * t371 + t68;
t267 = t268 + t400;
t255 = t261 ^ 2;
t250 = t255 * t266;
t236 = -0.2e1 * t245 + t316;
t233 = t249 + 0.2e1 * t309;
t214 = 0.2e1 * t321;
t157 = t203 * t258 - t204;
t99 = -t121 + t332;
t72 = t214 - t381 + t454;
t71 = t214 - t392 + (t274 + t310) * qJ(4) + t277;
t67 = -qJ(4) * t165 + t70;
t66 = (-t165 - t251) * pkin(3) + t283;
t53 = t258 * t87 + t259 * t86;
t41 = t258 * t70 + t259 * t69;
t40 = t258 * t69 - t259 * t70;
t39 = qJ(4) * t391 - t361 - t471;
t38 = t267 + t367;
t37 = qJ(4) * t100 - t362 - t433;
t30 = t370 * t391 + t362 - t470;
t29 = (t121 - t99) * pkin(5) + t267;
t28 = t268 + t367 + 0.2e1 * t400;
t27 = t370 * t100 - t361 - t434;
t22 = -qJ(6) * t111 + t26;
t21 = (-t111 - t374) * pkin(5) + t292;
t20 = -t260 * t29 + t306 * t99 - t433;
t19 = t471 + t263 * t28 + (-qJ(4) - t369) * t391;
t16 = -t263 * t29 + t299 * t99 - t434;
t15 = t470 - t260 * t28 + (-t368 - t370) * t391;
t14 = -pkin(8) * t17 - qJ(4) * t68;
t13 = t24 * t263 + t26 * t260;
t12 = t24 * t260 - t26 * t263;
t11 = -t17 + t448;
t10 = -pkin(8) * t18 - t370 * t68;
t9 = -t18 + t449;
t8 = -t21 * t260 + t22 * t263 + t448;
t7 = -t263 * t21 - t260 * t22 + t449;
t6 = t17 * t258 + t18 * t259;
t5 = -t17 * t259 + t18 * t258;
t4 = t12 * t258 + t13 * t259;
t3 = -t12 * t259 + t13 * t258;
t2 = -pkin(8) * t12 + (-t306 - t369) * t38;
t1 = -pkin(8) * t13 + (-t299 - t368) * t38;
t23 = [0, 0, 0, 0, 0, qJDD(1), t307, t296, 0, 0, t290 * t261, t233 * t264 + t236 * t261, t328 + t264 * (-t250 + t373), -t264 * t245 + t220, t261 * (t251 - t373) + t325, 0, t264 * t217 + pkin(1) * t236 + pkin(7) * (t264 * (-t251 - t373) - t328), -t261 * t217 - pkin(1) * t233 + pkin(7) * (-t325 - t261 * (-t250 - t373)), pkin(1) * (t250 + t251) + (t255 + t256) * t360 + t304, pkin(1) * t217 + pkin(7) * t304, t295, -t415, t414, t270, t446, t273, t261 * (t351 - t432) + t264 * (t108 - t435) + t410, t261 * (t350 - t463) + t264 * (t109 - t464) + t469, t288 * t117 + t261 * t289 + t445, pkin(7) * (t155 * t261 + t264 * t64) - t288 * t289, t295, t414, t415, t273, -t446, t270, t261 * (qJ(4) * t346 - t258 * t72 - t432) + t264 * (pkin(3) * t179 - qJ(4) * t383 - t435 + t70) + t410, t261 * (-qJ(3) * t117 - t258 * t66 + t259 * t67) + t264 * (-pkin(2) * t117 - pkin(3) * t275 - qJ(4) * t174) - pkin(1) * t117 + t445, t261 * (-pkin(3) * t342 + t259 * t71 + t463) + t264 * (t464 + qJ(4) * t279 + 0.2e1 * t322 + (t389 + t251) * pkin(3) - t387) - t469, t261 * (-qJ(3) * t40 + (pkin(3) * t258 - qJ(4) * t259) * t77) + t264 * (-pkin(2) * t40 + pkin(3) * t70 - qJ(4) * t69) - pkin(1) * t40 + pkin(7) * (t261 * t77 + t264 * t41), t297, t457, t438, t378, t472, t377, t261 * (-t258 * t27 + t259 * t37 - t453) + t264 * (-t31 + t450) - t455 + pkin(7) * (-t100 * t261 + t451), t261 * (-t258 * t30 + t259 * t39 - t476) + t264 * (-t32 - t475) - t478, t261 * (t11 * t259 - t258 * t9 - t461) + t264 * t460 + t456, t261 * (-qJ(3) * t5 - t10 * t258 + t14 * t259) + t264 * (-pkin(2) * t5 - qJ(4) * t18 + t370 * t17) - pkin(1) * t5 + pkin(7) * (t261 * t68 + t264 * t6), t297, t438, -t457, t377, -t472, t378, t261 * (-t16 * t258 + t20 * t259 - t453) + t264 * (t271 + t450) - t455 + pkin(7) * (-t261 * t99 + t451), t261 * (-t258 * t7 + t259 * t8 - t461) + t264 * (t363 + t460) + t456, t261 * (-t15 * t258 + t19 * t259 + t476) + t264 * (t276 + t475) + t478, t261 * (-qJ(3) * t3 - t1 * t258 + t2 * t259) + t264 * (-pkin(2) * t3 - qJ(4) * t13 + t370 * t12 + t364) - pkin(1) * t3 + pkin(7) * (t261 * t38 + t264 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t242, t250 - t251, t249, t242, t316, qJDD(2), -t198, -t199, 0, 0, t157, t382, t424, t285, t417, t284, -t350 + t421, -pkin(2) * t274 + t351 + t462, t447 + t64, -pkin(2) * t155 + qJ(3) * t64, t157, t424, -t382, t284, -t417, t285, qJ(4) * t347 + t259 * t72 + t421, t258 * t67 + t259 * t66 + t447, t258 * t71 + t274 * t311 - t462, qJ(3) * t41 + (-qJ(4) * t258 - t311) * t77, t53, -t444, t439, t380, -t467, t379, pkin(2) * t100 + t258 * t37 + t259 * t27 + t452, t258 * t39 + t259 * t30 - t474, t11 * t258 + t259 * t9 + t459, -pkin(2) * t68 + qJ(3) * t6 + t10 * t259 + t14 * t258, t53, t439, t444, t379, t467, t380, pkin(2) * t99 + t16 * t259 + t20 * t258 + t452, t258 * t8 + t259 * t7 + t459, t15 * t259 + t19 * t258 + t474, -pkin(2) * t38 + qJ(3) * t4 + t1 * t259 + t2 * t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t399, t274, t165, t155, 0, 0, 0, 0, 0, 0, -t399, t165, -t274, t77, 0, 0, 0, 0, 0, 0, -t100, -t391, -t111, t68, 0, 0, 0, 0, 0, 0, -t99, -t111, t391, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t275, t389, t70, 0, 0, 0, 0, 0, 0, t395, t78, t420, t17, 0, 0, 0, 0, 0, 0, t395, t420, -t78, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t336, t138, t385, -t336, -t101, t230, -t31, -t32, 0, 0, t336, t385, -t138, t230, t101, -t336, t271, t363, t276, t364; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t386, t385, t390, t26;];
tauJ_reg  = t23;
