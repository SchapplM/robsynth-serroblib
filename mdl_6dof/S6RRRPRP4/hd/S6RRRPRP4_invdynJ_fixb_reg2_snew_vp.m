% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRPRP4
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
% Datum: 2019-05-07 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRPRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:47:30
% EndTime: 2019-05-07 07:48:01
% DurationCPUTime: 9.70s
% Computational Cost: add. (22718->478), mult. (47246->567), div. (0->0), fcn. (32825->8), ass. (0->300)
t298 = cos(qJ(3));
t290 = qJD(2) + qJD(3);
t288 = t290 ^ 2;
t294 = sin(qJ(3));
t299 = cos(qJ(2));
t295 = sin(qJ(2));
t369 = t295 * t298;
t265 = (t294 * t299 + t369) * qJD(1);
t413 = t265 ^ 2;
t425 = -t413 - t288;
t289 = qJDD(2) + qJDD(3);
t364 = qJD(1) * t295;
t263 = -t298 * t299 * qJD(1) + t294 * t364;
t375 = t265 * t263;
t461 = t375 + t289;
t485 = t461 * t294;
t176 = -t298 * t425 + t485;
t484 = t461 * t298;
t178 = t294 * t425 + t484;
t515 = pkin(7) * (t176 * t295 - t178 * t299);
t514 = pkin(2) * t176;
t513 = pkin(8) * t176;
t512 = pkin(8) * t178;
t293 = sin(qJ(5));
t297 = cos(qJ(5));
t239 = t263 * t293 + t290 * t297;
t283 = t295 * qJDD(1);
t357 = qJD(1) * qJD(2);
t350 = t299 * t357;
t272 = t283 + t350;
t284 = t299 * qJDD(1);
t351 = t295 * t357;
t273 = t284 - t351;
t343 = t294 * t272 - t298 * t273;
t203 = qJD(3) * t265 + t343;
t346 = -t297 * t203 + t293 * t289;
t258 = qJD(5) + t265;
t359 = qJD(5) - t258;
t323 = -t239 * t359 - t346;
t415 = t258 ^ 2;
t237 = -t297 * t263 + t290 * t293;
t416 = t237 ^ 2;
t211 = t416 - t415;
t335 = t298 * t272 + t294 * t273;
t204 = -qJD(3) * t263 + t335;
t319 = -qJDD(5) - t204;
t380 = t239 * t237;
t143 = -t380 + t319;
t392 = t143 * t297;
t99 = t211 * t293 - t392;
t510 = t295 * (t294 * t99 + t323 * t298) - t299 * (-t323 * t294 + t298 * t99);
t235 = t239 ^ 2;
t428 = -t235 - t415;
t95 = -t293 * t428 + t392;
t509 = pkin(1) * t95;
t508 = pkin(2) * t95;
t393 = t143 * t293;
t93 = t297 * t428 + t393;
t507 = t294 * t93;
t506 = t298 * t93;
t411 = pkin(3) + pkin(9);
t505 = t411 * t93;
t504 = t411 * t95;
t501 = -pkin(4) * t93 + qJ(4) * t95;
t247 = -t413 + t288;
t219 = t375 - t289;
t497 = t219 * t298;
t498 = t219 * t294;
t500 = t295 * (t247 * t294 + t497) + t299 * (-t247 * t298 + t498);
t414 = t263 ^ 2;
t246 = t414 - t288;
t499 = t295 * (-t246 * t298 + t485) - t299 * (t246 * t294 + t484);
t358 = qJD(5) + t258;
t322 = t239 * t358 + t346;
t426 = -t380 - t319;
t391 = t426 * t297;
t422 = -t415 - t416;
t435 = t293 * t422 + t391;
t454 = t294 * t322 - t298 * t435;
t492 = pkin(2) * t454;
t155 = -t416 - t235;
t336 = t293 * t203 + t297 * t289;
t130 = t237 * t359 - t336;
t419 = t130 * t297 + t293 * t323;
t457 = t155 * t294 - t298 * t419;
t491 = pkin(2) * t457;
t218 = -t288 - t414;
t157 = t218 * t294 - t497;
t160 = -t218 * t298 - t498;
t490 = pkin(7) * (t157 * t295 + t160 * t299);
t489 = pkin(8) * t454;
t488 = pkin(8) * t457;
t136 = t293 * t426;
t212 = -t235 + t415;
t462 = -t212 * t297 - t136;
t487 = t294 * t462;
t486 = t298 * t462;
t103 = -t211 * t297 - t393;
t420 = -t293 * t130 + t297 * t323;
t456 = t155 * t298 + t294 * t419;
t483 = -pkin(2) * t420 + pkin(8) * t456;
t436 = t297 * t422 - t136;
t455 = t294 * t435 + t298 * t322;
t482 = -pkin(2) * t436 + pkin(8) * t455;
t481 = pkin(7) * (-t295 * t454 + t299 * t455) - pkin(1) * t436;
t480 = pkin(7) * (-t295 * t457 + t299 * t456) - pkin(1) * t420;
t478 = pkin(2) * t157;
t476 = pkin(8) * t157;
t475 = pkin(8) * t160;
t151 = -qJD(5) * t237 + t336;
t214 = t237 * t258;
t128 = t151 - t214;
t474 = qJ(6) * t128;
t361 = qJD(3) - t290;
t175 = t263 * t361 - t335;
t473 = t175 * t294;
t472 = t175 * t298;
t186 = t235 - t416;
t471 = t186 * t294;
t470 = t186 * t298;
t463 = t411 * t436;
t460 = -t212 * t293 + t391;
t459 = pkin(4) * t155 - t411 * t420;
t458 = qJ(4) * t155 - t411 * t419;
t453 = qJ(4) * t322 - t411 * t435;
t452 = pkin(4) * t435 - qJ(4) * t436;
t450 = 2 * qJD(4);
t423 = -t414 - t413;
t448 = pkin(1) * t423;
t446 = pkin(2) * t423;
t292 = t299 ^ 2;
t301 = qJD(1) ^ 2;
t296 = sin(qJ(1));
t410 = cos(qJ(1));
t347 = g(1) * t296 - t410 * g(2);
t330 = qJDD(1) * pkin(1) + t347;
t331 = qJD(2) * pkin(2) - pkin(8) * t364;
t207 = pkin(2) * t273 + (pkin(8) * t292 + pkin(7)) * t301 - t331 * t364 + t330;
t373 = t290 * t265;
t438 = pkin(3) * t373 - t265 * t450 - t207;
t42 = pkin(4) * t419 - qJ(4) * t420;
t376 = t263 * t290;
t368 = t295 * t301;
t332 = g(1) * t410 + t296 * g(2);
t400 = qJDD(1) * pkin(7);
t268 = -t301 * pkin(1) - t332 + t400;
t374 = t268 * t295;
t185 = qJDD(2) * pkin(2) - pkin(8) * t272 - t374 + (pkin(2) * t368 + pkin(8) * t357 - g(3)) * t299;
t244 = -t295 * g(3) + t299 * t268;
t286 = t292 * t301;
t188 = -pkin(2) * t286 + t273 * pkin(8) - qJD(2) * t331 + t244;
t134 = -t298 * t185 + t188 * t294;
t225 = pkin(3) * t263 - qJ(4) * t265;
t90 = -t289 * pkin(3) - t288 * qJ(4) + t265 * t225 + qJDD(4) + t134;
t305 = t219 * pkin(9) + (t204 + t376) * pkin(4) + t90;
t245 = pkin(4) * t265 - pkin(9) * t290;
t345 = -t204 + t376;
t433 = qJ(4) * t345;
t302 = t433 + t438;
t54 = -pkin(4) * t414 + t203 * t411 - t245 * t265 + t302;
t37 = t293 * t54 - t297 * t305;
t38 = t293 * t305 + t297 * t54;
t19 = t293 * t38 - t297 * t37;
t360 = qJD(3) + t290;
t316 = -t263 * t360 + t335;
t432 = t294 * t316;
t430 = t298 * t316;
t424 = t413 - t414;
t150 = -qJD(5) * t239 - t346;
t421 = -t150 * pkin(5) - t474;
t377 = t258 * t297;
t209 = t239 * t377;
t378 = t258 * t293;
t354 = t237 * t378;
t340 = t209 + t354;
t418 = t295 * (-t294 * t340 - t298 * t319) + t299 * (-t294 * t319 + t298 * t340);
t327 = -t150 * t297 - t354;
t355 = t298 * t380;
t356 = t294 * t380;
t417 = t295 * (-t294 * t327 - t355) + t299 * (t298 * t327 - t356);
t412 = 2 * qJD(6);
t409 = pkin(3) * t294;
t408 = pkin(3) * t298;
t407 = pkin(5) * t297;
t184 = pkin(5) * t237 - qJ(6) * t239;
t339 = -qJ(6) * t319 - t237 * t184 + t258 * t412 + t38;
t29 = -pkin(5) * t415 + t339;
t31 = pkin(5) * t319 - qJ(6) * t415 + t184 * t239 + qJDD(6) + t37;
t406 = -pkin(5) * t31 + qJ(6) * t29;
t135 = t185 * t294 + t188 * t298;
t312 = -t288 * pkin(3) + qJ(4) * t289 - t225 * t263 + t135;
t362 = qJD(4) * t290;
t82 = t312 + 0.2e1 * t362;
t405 = -pkin(3) * t90 + qJ(4) * t82;
t311 = -t203 * pkin(4) - pkin(9) * t414 + t312;
t63 = (t450 + t245) * t290 + t311;
t403 = t293 * t63;
t76 = -t134 * t298 + t135 * t294;
t402 = t295 * t76;
t61 = t297 * t63;
t401 = qJ(6) * t297;
t399 = t322 * t293;
t398 = t322 * t297;
t397 = t128 * t293;
t390 = t207 * t294;
t389 = t207 * t298;
t379 = t239 * t258;
t372 = t290 * t294;
t371 = t290 * t298;
t278 = t299 * t368;
t370 = t295 * (qJDD(2) + t278);
t367 = t299 * (qJDD(2) - t278);
t366 = pkin(5) * t130 + qJ(6) * t323;
t168 = -t265 * t361 - t343;
t365 = pkin(3) * t175 + qJ(4) * t168;
t353 = t237 * t377;
t349 = qJ(4) * t294 + pkin(2);
t348 = qJ(6) * t293 + pkin(4);
t77 = t134 * t294 + t298 * t135;
t243 = g(3) * t299 + t374;
t344 = t295 * t243 + t299 * t244;
t342 = qJ(4) * t63 - t19 * t411;
t208 = t239 * t378;
t341 = t208 - t353;
t20 = t293 * t37 + t297 * t38;
t317 = -t237 * t358 + t336;
t337 = qJ(4) * t317 - t505 + t61;
t303 = t239 * t412 - t245 * t290 - t311 - 0.2e1 * t362 - t421;
t34 = -pkin(5) * t379 + t303 + t474;
t329 = -pkin(5) * t397 - qJ(4) * t128 + t297 * t34 + t505;
t328 = -t150 * t293 + t353;
t326 = t403 + t453;
t325 = -t19 + t458;
t23 = (-t155 - t415) * pkin(5) + t339;
t25 = -qJ(6) * t155 + t31;
t324 = -t23 * t293 + t297 * t25 + t458;
t10 = t29 * t293 - t297 * t31;
t41 = (pkin(5) * t258 - (2 * qJD(6))) * t239 + t63 + t421;
t321 = -t411 * t10 + (pkin(5) * t293 + qJ(4) - t401) * t41;
t320 = -pkin(5) * t428 - qJ(6) * t143 + t29;
t315 = pkin(3) * t219 - qJ(4) * t218 + t90;
t35 = t303 + (-t322 - t379) * pkin(5);
t314 = -t293 * t35 - t322 * t401 + t453;
t313 = pkin(5) * t426 + qJ(6) * t422 - t31;
t310 = t295 * (t298 * t204 - t265 * t372) + t299 * (t294 * t204 + t265 * t371);
t114 = -t151 * t293 - t209;
t309 = t295 * (-t114 * t294 + t355) + t299 * (t298 * t114 + t356);
t308 = t295 * (t203 * t294 + t263 * t371) + t299 * (-t298 * t203 + t263 * t372);
t307 = -pkin(3) * t425 + qJ(4) * t461 + t82;
t306 = (t295 * (-t263 * t298 + t265 * t294) + t299 * (-t263 * t294 - t265 * t298)) * t290;
t304 = -pkin(3) * t203 - t438;
t300 = qJD(2) ^ 2;
t291 = t295 ^ 2;
t285 = t291 * t301;
t274 = t284 - 0.2e1 * t351;
t271 = t283 + 0.2e1 * t350;
t267 = pkin(7) * t301 + t330;
t169 = -t203 + t373;
t167 = t203 + t373;
t166 = t265 * t360 + t343;
t132 = t151 + t214;
t115 = t151 * t297 - t208;
t108 = t169 * t298 - t473;
t107 = t168 * t298 - t473;
t106 = t169 * t294 + t472;
t105 = t168 * t294 + t472;
t79 = -qJ(4) * t423 + t90;
t78 = -pkin(3) * t423 + t82;
t73 = -t293 * t317 - t398;
t72 = t397 + t398;
t69 = -t297 * t317 + t399;
t68 = t128 * t297 - t399;
t65 = (t167 + t203) * pkin(3) + t302;
t64 = (t316 - t345) * qJ(4) + t304;
t58 = t298 * t317 + t507;
t56 = t294 * t317 - t506;
t52 = -t128 * t298 - t507;
t50 = -t128 * t294 + t506;
t43 = t294 * t82 - t298 * t90;
t32 = pkin(4) * t317 - t403 - t504;
t28 = pkin(4) * t322 - t463 + t61;
t26 = t366 + t42;
t22 = -t38 - t501;
t21 = -t37 + t452;
t17 = t313 + t452;
t16 = t320 + t501;
t15 = -t297 * t35 + t322 * t348 - t463;
t14 = -t293 * t34 + t504 + (-pkin(4) - t407) * t128;
t13 = t19 * t294 + t298 * t63;
t12 = -t19 * t298 + t294 * t63;
t11 = t29 * t297 + t293 * t31;
t8 = -t20 + t459;
t7 = t10 * t294 + t298 * t41;
t6 = -t10 * t298 + t294 * t41;
t5 = -t23 * t297 - t25 * t293 + t459;
t4 = pkin(4) * t19 - qJ(4) * t20;
t3 = pkin(4) * t63 - t20 * t411;
t2 = pkin(4) * t10 - qJ(4) * t11 + t406;
t1 = -t411 * t11 + (t348 + t407) * t41;
t9 = [0, 0, 0, 0, 0, qJDD(1), t347, t332, 0, 0, (t272 + t350) * t295, t271 * t299 + t274 * t295, t370 + t299 * (-t285 + t300), (t273 - t351) * t299, t295 * (t286 - t300) + t367, 0, t299 * t267 + pkin(1) * t274 + pkin(7) * (t299 * (-t286 - t300) - t370), -t295 * t267 - pkin(1) * t271 + pkin(7) * (-t367 - t295 * (-t285 - t300)), pkin(1) * (t285 + t286) + (t291 + t292) * t400 + t344, pkin(1) * t267 + pkin(7) * t344, t310, t295 * (-t166 * t298 - t432) + t299 * (-t166 * t294 + t430), -t500, t308, -t499, t306, t295 * (-t390 - t476) + t299 * (-pkin(2) * t166 + t389 - t475) - pkin(1) * t166 - t490, t295 * (-t389 + t513) + t299 * (-pkin(2) * t316 - t390 - t512) - pkin(1) * t316 + t515, t295 * (-pkin(8) * t106 - t76) + t299 * (pkin(8) * t108 - t446 + t77) - t448 + pkin(7) * (-t106 * t295 + t108 * t299), -pkin(8) * t402 + t299 * (pkin(2) * t207 + pkin(8) * t77) + pkin(1) * t207 + pkin(7) * (t299 * t77 - t402), t306, t500, t499, t310, t295 * (-t167 * t298 - t432) + t299 * (-t167 * t294 + t430), t308, t295 * (-pkin(8) * t105 - t294 * t78 + t298 * t79) + t299 * (pkin(8) * t107 + t294 * t79 + t298 * t78 - t446) - t448 + pkin(7) * (-t105 * t295 + t107 * t299), t295 * (-t294 * t65 + t476) + t299 * (t298 * t65 + t475) + t490 + (qJ(4) * t369 + t299 * t349 + pkin(1)) * t167, t295 * (t298 * t64 - t513) + t299 * (t294 * t64 + t512) - t515 + (-t295 * t409 + t299 * (pkin(2) + t408) + pkin(1)) * t316, (t295 * (qJ(4) * t298 - t409) + t299 * (t349 + t408) + pkin(1)) * (t304 - t433) + (pkin(7) + pkin(8)) * (-t295 * t43 + t299 * (t294 * t90 + t298 * t82)), t309, t295 * (-t294 * t69 + t470) + t299 * (t298 * t69 + t471), t295 * (-t130 * t298 - t487) + t299 * (-t130 * t294 + t486), t417, t510, t418, t295 * (t21 * t298 - t28 * t294 - t489) + t299 * (t21 * t294 + t298 * t28 + t482) + t481, t295 * (-pkin(8) * t56 + t22 * t298 - t294 * t32) + t299 * (pkin(8) * t58 + t22 * t294 + t298 * t32 - t508) - t509 + pkin(7) * (-t295 * t56 + t299 * t58), t295 * (-t294 * t8 + t298 * t42 - t488) + t299 * (t294 * t42 + t298 * t8 + t483) + t480, t295 * (-pkin(8) * t12 - t294 * t3 + t298 * t4) + t299 * (-pkin(2) * t20 + pkin(8) * t13 + t294 * t4 + t298 * t3) - pkin(1) * t20 + pkin(7) * (-t12 * t295 + t13 * t299), t309, t295 * (t132 * t298 - t487) + t299 * (t132 * t294 + t486), t295 * (-t294 * t68 - t470) + t299 * (t298 * t68 - t471), t418, -t510, t417, t295 * (-t15 * t294 + t17 * t298 - t489) + t299 * (t15 * t298 + t17 * t294 + t482) + t481, t295 * (t26 * t298 - t294 * t5 - t488) + t299 * (t26 * t294 + t298 * t5 + t483) + t480, t295 * (-pkin(8) * t50 - t14 * t294 + t16 * t298) + t299 * (pkin(8) * t52 + t14 * t298 + t16 * t294 + t508) + t509 + pkin(7) * (-t295 * t50 + t299 * t52), t295 * (-pkin(8) * t6 - t1 * t294 + t2 * t298) + t299 * (-pkin(2) * t11 + pkin(8) * t7 + t1 * t298 + t2 * t294) - pkin(1) * t11 + pkin(7) * (-t295 * t6 + t299 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, t285 - t286, t283, t278, t284, qJDD(2), -t243, -t244, 0, 0, t375, t424, -t175, -t375, t168, t289, -t134 + t478, -t135 - t514, pkin(2) * t106, pkin(2) * t76, t289, t175, -t169, t375, t424, -t375, pkin(2) * t105 + t365, t315 - t478, t307 + t514, pkin(2) * t43 + t405, t115, t73, t460, t328, -t103, t341, t326 + t492, pkin(2) * t56 + t337, t325 + t491, pkin(2) * t12 + t342, t115, t460, t72, t341, t103, t328, t314 + t492, t324 + t491, pkin(2) * t50 + t329, pkin(2) * t6 + t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t375, t424, -t175, -t375, t168, t289, -t134, -t135, 0, 0, t289, t175, -t169, t375, t424, -t375, t365, t315, t307, t405, t115, t73, t460, t328, -t103, t341, t326, t337, t325, t342, t115, t460, t72, t341, t103, t328, t314, t324, t329, t321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, -t219, t425, t90, 0, 0, 0, 0, 0, 0, t435, t93, t419, t19, 0, 0, 0, 0, 0, 0, t435, t419, -t93, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t380, t186, -t130, -t380, t323, -t319, -t37, -t38, 0, 0, t380, t132, -t186, -t319, -t323, -t380, t313, t366, t320, t406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, -t130, t428, t31;];
tauJ_reg  = t9;
