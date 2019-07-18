% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:32:23
% EndTime: 2019-05-05 21:32:51
% DurationCPUTime: 15.64s
% Computational Cost: add. (45288->533), mult. (109496->717), div. (0->0), fcn. (84018->10), ass. (0->343)
t315 = sin(pkin(9));
t317 = cos(pkin(9));
t320 = sin(qJ(3));
t322 = cos(qJ(3));
t366 = t317 * qJDD(1);
t367 = t315 * qJDD(1);
t296 = t320 * t366 + t322 * t367;
t374 = qJD(1) * t317;
t382 = t315 * t320;
t297 = qJD(1) * t382 - t322 * t374;
t369 = t297 * qJD(3);
t272 = -t369 + t296;
t381 = t317 * t320;
t299 = (t315 * t322 + t381) * qJD(1);
t319 = sin(qJ(4));
t321 = cos(qJ(4));
t282 = qJD(3) * t319 + t299 * t321;
t232 = -t282 * qJD(4) + t321 * qJDD(3) - t319 * t272;
t281 = -t321 * qJD(3) + t299 * t319;
t344 = -t319 * qJDD(3) - t321 * t272;
t233 = -qJD(4) * t281 - t344;
t314 = sin(pkin(10));
t316 = cos(pkin(10));
t345 = t232 * t314 + t233 * t316;
t245 = t316 * t281 + t282 * t314;
t292 = qJD(4) + t297;
t393 = t245 * t292;
t441 = t345 - t393;
t347 = t320 * t367 - t322 * t366;
t368 = t299 * qJD(3);
t271 = -t347 - t368;
t262 = qJDD(4) - t271;
t247 = -t281 * t314 + t282 * t316;
t392 = t247 * t245;
t181 = -t392 - t262;
t399 = t181 * t314;
t244 = t247 ^ 2;
t425 = t292 ^ 2;
t443 = -t244 - t425;
t128 = t316 * t443 + t399;
t398 = t181 * t316;
t130 = -t314 * t443 + t398;
t78 = t128 * t319 - t130 * t321;
t62 = t320 * t78 + t322 * t441;
t64 = -t320 * t441 + t322 * t78;
t81 = t128 * t321 + t130 * t319;
t525 = pkin(1) * t81 - qJ(2) * (t315 * t62 - t317 * t64);
t524 = pkin(7) * t62;
t523 = pkin(2) * t81 + pkin(7) * t64;
t184 = t316 * t232 - t233 * t314;
t387 = t292 * t247;
t339 = t184 + t387;
t427 = t245 ^ 2;
t220 = t427 - t425;
t134 = -t220 * t314 + t398;
t138 = -t220 * t316 - t399;
t91 = t134 * t319 - t138 * t321;
t521 = t315 * (t320 * t339 + t322 * t91) + t317 * (t320 * t91 - t339 * t322);
t518 = pkin(3) * t81;
t517 = pkin(8) * t81;
t516 = pkin(3) * t441 + pkin(8) * t78;
t195 = t427 - t244;
t338 = -t184 + t387;
t103 = -t338 * t314 + t316 * t441;
t404 = t441 * t314;
t106 = t338 * t316 + t404;
t55 = t103 * t319 + t106 * t321;
t513 = t315 * (t320 * t195 + t322 * t55) + t317 * (-t195 * t322 + t320 * t55);
t510 = pkin(4) * t128;
t509 = qJ(5) * t128;
t508 = qJ(5) * t130;
t507 = t134 * t321 + t138 * t319;
t440 = -t392 + t262;
t397 = t440 * t314;
t438 = -t425 - t427;
t454 = t316 * t438 - t397;
t171 = t316 * t440;
t455 = t314 * t438 + t171;
t471 = -t319 * t455 + t321 * t454;
t492 = t320 * t471 - t322 * t338;
t506 = pkin(7) * t492;
t470 = t319 * t454 + t321 * t455;
t491 = t320 * t338 + t322 * t471;
t505 = -pkin(2) * t470 + pkin(7) * t491;
t504 = -t103 * t321 + t106 * t319;
t503 = qJ(2) * (-t315 * t492 + t317 * t491) - pkin(1) * t470;
t439 = t393 + t345;
t221 = -t244 + t425;
t476 = t316 * t221 + t397;
t477 = -t221 * t314 + t171;
t489 = -t319 * t476 + t321 * t477;
t502 = t315 * (t320 * t439 + t322 * t489) + t317 * (t320 * t489 - t322 * t439);
t499 = pkin(3) * t470;
t162 = -t427 - t244;
t433 = t314 * t439 + t316 * t339;
t434 = t314 * t339 - t316 * t439;
t452 = -t319 * t434 + t321 * t433;
t473 = -t162 * t322 + t320 * t452;
t498 = pkin(7) * t473;
t497 = pkin(8) * t470;
t451 = t319 * t433 + t321 * t434;
t472 = t320 * t162 + t322 * t452;
t494 = -pkin(2) * t451 + pkin(7) * t472;
t493 = -pkin(3) * t338 + pkin(8) * t471;
t490 = t319 * t477 + t321 * t476;
t488 = qJ(2) * (-t315 * t473 + t317 * t472) - pkin(1) * t451;
t485 = pkin(3) * t451;
t417 = pkin(4) * t455;
t484 = pkin(8) * t451;
t483 = qJ(5) * t454;
t482 = qJ(5) * t455;
t475 = -pkin(4) * t162 + qJ(5) * t433;
t474 = -pkin(3) * t162 + pkin(8) * t452;
t418 = pkin(4) * t434;
t467 = qJ(5) * t434;
t466 = qJ(6) * t441;
t251 = t282 * t281;
t442 = -t251 + t262;
t463 = t319 * t442;
t275 = t299 * t297;
t437 = -t275 + qJDD(3);
t462 = t320 * t437;
t459 = t321 * t442;
t458 = t322 * t437;
t419 = sin(qJ(1));
t420 = cos(qJ(1));
t350 = g(1) * t419 - t420 * g(2);
t342 = -qJDD(2) + t350;
t424 = qJD(1) ^ 2;
t373 = t424 * qJ(2);
t410 = qJDD(1) * pkin(1);
t293 = t342 + t373 + t410;
t311 = t315 ^ 2;
t312 = t317 ^ 2;
t375 = t311 + t312;
t453 = t373 * t375 - t293 - t410;
t360 = pkin(2) * t317 + pkin(1);
t266 = t360 * qJDD(1) + (pkin(7) * t375 + qJ(2)) * t424 + t342;
t362 = t320 * t392;
t365 = t322 * t392;
t385 = t292 * t316;
t364 = t245 * t385;
t340 = -t184 * t314 + t364;
t386 = t292 * t314;
t349 = t316 * t184 + t245 * t386;
t429 = -t319 * t349 + t321 * t340;
t450 = t315 * (t322 * t429 - t362) + t317 * (t320 * t429 + t365);
t258 = t322 * t262;
t379 = t320 * t262;
t336 = (-t245 * t314 - t247 * t316) * t292;
t216 = t247 * t386;
t348 = t216 - t364;
t430 = -t319 * t336 + t321 * t348;
t449 = t315 * (t322 * t430 + t379) + t317 * (t320 * t430 - t258);
t447 = pkin(7) + qJ(2);
t257 = t292 * t281;
t204 = t233 + t257;
t372 = qJD(5) * t245;
t239 = -0.2e1 * t372;
t370 = qJD(6) * t292;
t444 = t239 + 0.2e1 * t370;
t194 = pkin(5) * t245 - qJ(6) * t247;
t263 = pkin(3) * t297 - pkin(8) * t299;
t334 = g(1) * t420 + g(2) * t419;
t326 = (-t447 * qJDD(1) + (qJD(1) * t360 - (2 * qJD(2))) * qJD(1) + t334) * t315;
t415 = t317 * g(3);
t325 = t326 - t415;
t331 = qJDD(1) * qJ(2) - t334;
t422 = 2 * qJD(2);
t351 = -g(3) * t315 + t317 * (-pkin(1) * t424 + t331) + t374 * t422;
t256 = -pkin(2) * t312 * t424 + pkin(7) * t366 + t351;
t377 = t322 * t256;
t423 = qJD(3) ^ 2;
t175 = -pkin(3) * t423 + qJDD(3) * pkin(8) - t297 * t263 + t320 * t325 + t377;
t187 = (-t272 + t369) * pkin(8) + (-t271 + t368) * pkin(3) - t266;
t118 = t319 * t175 - t321 * t187;
t87 = pkin(4) * t442 - qJ(5) * t204 - t118;
t119 = t321 * t175 + t319 * t187;
t253 = pkin(4) * t292 - qJ(5) * t282;
t426 = t281 ^ 2;
t96 = -pkin(4) * t426 + qJ(5) * t232 - t253 * t292 + t119;
t414 = t314 * t87 + t316 * t96;
t353 = t262 * qJ(6) - t245 * t194 + t414;
t436 = -t510 - pkin(5) * (t443 + t425) - qJ(6) * t181 + t353;
t431 = t319 * t348 + t321 * t336;
t428 = t319 * t340 + t321 * t349;
t280 = t282 ^ 2;
t294 = t297 ^ 2;
t295 = t299 ^ 2;
t357 = t314 * t96 - t316 * t87;
t371 = qJD(5) * t247;
t51 = t357 + 0.2e1 * t371;
t52 = t239 + t414;
t28 = t314 * t52 - t316 * t51;
t421 = pkin(4) * t28;
t416 = pkin(5) * t316;
t413 = t28 * t319;
t412 = t28 * t321;
t411 = qJ(6) * t316;
t207 = t320 * t256 - t322 * t325;
t174 = -qJDD(3) * pkin(3) - t423 * pkin(8) + t299 * t263 + t207;
t110 = -t232 * pkin(4) - t426 * qJ(5) + t282 * t253 + qJDD(5) + t174;
t409 = t110 * t314;
t408 = t110 * t316;
t208 = -g(3) * t381 + t320 * t326 + t377;
t159 = -t207 * t322 + t320 * t208;
t402 = t159 * t315;
t401 = t174 * t319;
t400 = t174 * t321;
t213 = t251 + t262;
t396 = t213 * t319;
t395 = t213 * t321;
t391 = t266 * t320;
t390 = t266 * t322;
t267 = qJDD(3) + t275;
t389 = t267 * t322;
t388 = t282 * t292;
t384 = t292 * t319;
t383 = t292 * t321;
t378 = t320 * t267;
t363 = t322 * t251;
t361 = t320 * t251;
t359 = -pkin(3) * t322 - pkin(2);
t358 = -qJ(6) * t314 - pkin(4);
t29 = t314 * t51 + t316 * t52;
t75 = t118 * t319 + t321 * t119;
t160 = t207 * t320 + t322 * t208;
t355 = t315 * (t415 + ((-pkin(1) * qJD(1) + t422) * qJD(1) + t331) * t315) + t317 * t351;
t354 = (0.2e1 * qJD(5) + t194) * t247;
t352 = -t414 + t510;
t74 = -t118 * t321 + t119 * t319;
t343 = t232 + t388;
t341 = t353 + t444;
t337 = -t262 * pkin(5) - qJ(6) * t425 + qJDD(6) + t357;
t42 = -pkin(5) * t425 + t341;
t43 = t354 + t337;
t22 = t314 * t42 - t316 * t43;
t335 = pkin(4) * t22 - pkin(5) * t43 + qJ(6) * t42;
t333 = -pkin(5) * t439 + qJ(6) * t339 + t418;
t147 = t247 * t385 + t314 * t345;
t148 = t316 * t345 - t216;
t100 = -t147 * t319 + t148 * t321;
t330 = t315 * (t322 * t100 + t362) + t317 * (t320 * t100 - t365);
t329 = -pkin(5) * t440 - qJ(6) * t438 + t337 - t417;
t328 = -t184 * pkin(5) + t110 - t466;
t327 = 0.2e1 * qJD(6) * t247 - t328;
t307 = t312 * qJDD(1);
t306 = t311 * qJDD(1);
t300 = t375 * t424;
t286 = -t295 - t423;
t285 = -t295 + t423;
t284 = t294 - t423;
t273 = 0.2e1 * t369 - t296;
t270 = t347 + 0.2e1 * t368;
t264 = -t423 - t294;
t255 = -t280 + t425;
t254 = -t425 + t426;
t250 = -t294 - t295;
t249 = t280 - t426;
t241 = -0.2e1 * t371;
t240 = 0.2e1 * t372;
t237 = -t280 - t425;
t236 = -t320 * t286 - t389;
t235 = t286 * t322 - t378;
t234 = -t425 - t426;
t225 = t280 + t426;
t224 = t320 * t296 - t322 * t347;
t223 = -t296 * t322 - t320 * t347;
t218 = t264 * t322 - t462;
t217 = t320 * t264 + t458;
t211 = (-t281 * t321 + t282 * t319) * t292;
t205 = (qJD(4) + t292) * t281 + t344;
t203 = t233 - t257;
t201 = t232 - t388;
t193 = t233 * t321 - t282 * t384;
t192 = -t232 * t319 + t281 * t383;
t189 = t254 * t321 - t396;
t188 = -t255 * t319 + t459;
t173 = -t237 * t319 - t395;
t172 = t237 * t321 - t396;
t165 = t234 * t321 - t463;
t164 = t234 * t319 + t459;
t142 = t204 * t319 + t321 * t343;
t141 = t201 * t321 - t203 * t319;
t140 = -t204 * t321 + t319 * t343;
t127 = t173 * t322 - t320 * t205;
t126 = t320 * t173 + t205 * t322;
t125 = t165 * t322 - t320 * t201;
t124 = t320 * t165 + t201 * t322;
t116 = -pkin(8) * t172 + t400;
t115 = t142 * t322 - t320 * t225;
t114 = t320 * t142 + t225 * t322;
t113 = -pkin(8) * t164 + t401;
t101 = -pkin(3) * t172 + t119;
t97 = t147 * t321 + t148 * t319;
t95 = -pkin(3) * t164 + t118;
t80 = t408 - t509;
t73 = t409 - t482;
t66 = -pkin(4) * t441 + t409 + t508;
t61 = (pkin(5) * t292 - 0.2e1 * qJD(6)) * t247 + t328;
t60 = -pkin(4) * t338 - t408 + t483;
t59 = -pkin(8) * t140 - t74;
t49 = t327 + (-t338 - t387) * pkin(5);
t48 = -pkin(5) * t387 + t327 + t466;
t41 = -t418 - t485;
t40 = -qJ(6) * t162 + t43;
t39 = (-t162 - t425) * pkin(5) + t341;
t38 = -t314 * t49 - t338 * t411 - t482;
t37 = -pkin(5) * t404 + t316 * t48 + t509;
t36 = t316 * t49 + t338 * t358 + t483;
t35 = -t508 + t314 * t48 + (pkin(4) + t416) * t441;
t34 = t239 - t352 - t518;
t33 = -t333 - t485;
t32 = -t319 * t66 + t321 * t80 - t517;
t31 = -t417 + t51 - t499;
t30 = -t319 * t60 + t321 * t73 - t497;
t27 = t329 + t354 - t499;
t26 = -pkin(4) * t110 + qJ(5) * t29;
t25 = t240 - 0.2e1 * t370 - t436 + t518;
t24 = -t28 - t467;
t23 = t314 * t43 + t316 * t42;
t21 = t29 + t475;
t20 = -t314 * t39 + t316 * t40 - t467;
t19 = t314 * t40 + t316 * t39 + t475;
t18 = -t319 * t36 + t321 * t38 - t497;
t17 = -t319 * t35 + t321 * t37 + t517;
t16 = t29 * t321 - t413;
t15 = t29 * t319 + t412;
t14 = -qJ(5) * t22 + (pkin(5) * t314 - t411) * t61;
t13 = t320 * t110 + t16 * t322;
t12 = -t110 * t322 + t320 * t16;
t11 = qJ(5) * t23 + (t358 - t416) * t61;
t10 = -t22 * t319 + t23 * t321;
t9 = t22 * t321 + t23 * t319;
t8 = t10 * t322 + t320 * t61;
t7 = t320 * t10 - t322 * t61;
t6 = -t21 * t319 + t24 * t321 - t484;
t5 = -pkin(3) * t15 - t421;
t4 = -t19 * t319 + t20 * t321 - t484;
t3 = -pkin(8) * t15 - qJ(5) * t412 - t26 * t319;
t2 = -pkin(3) * t9 - t335;
t1 = -pkin(8) * t9 - t11 * t319 + t14 * t321;
t44 = [0, 0, 0, 0, 0, qJDD(1), t350, t334, 0, 0, t306, 0.2e1 * t315 * t366, 0, t307, 0, 0, -t453 * t317, t453 * t315, pkin(1) * t300 + qJ(2) * (t307 + t306) + t355, pkin(1) * t293 + qJ(2) * t355, t315 * (t272 * t322 - t320 * t368) + t317 * (t320 * t272 + t322 * t368), t315 * (-t270 * t322 + t320 * t273) + t317 * (-t320 * t270 - t273 * t322), t315 * (-t320 * t285 + t458) + t317 * (t285 * t322 + t462), t315 * (-t320 * t271 + t322 * t369) + t317 * (t271 * t322 + t320 * t369), t315 * (t284 * t322 - t378) + t317 * (t320 * t284 + t389), (t315 * (-t297 * t322 + t299 * t320) + t317 * (-t297 * t320 - t299 * t322)) * qJD(3), t315 * (-pkin(7) * t217 - t391) + t317 * (-pkin(2) * t270 + pkin(7) * t218 + t390) - pkin(1) * t270 + qJ(2) * (-t217 * t315 + t218 * t317), t315 * (-pkin(7) * t235 - t390) + t317 * (pkin(2) * t273 + pkin(7) * t236 - t391) + pkin(1) * t273 + qJ(2) * (-t235 * t315 + t236 * t317), t315 * (-pkin(7) * t223 - t159) + t317 * (-pkin(2) * t250 + pkin(7) * t224 + t160) - pkin(1) * t250 + qJ(2) * (-t223 * t315 + t224 * t317), -pkin(7) * t402 + t317 * (pkin(2) * t266 + pkin(7) * t160) + pkin(1) * t266 + qJ(2) * (t160 * t317 - t402), t315 * (t193 * t322 + t361) + t317 * (t320 * t193 - t363), t315 * (t141 * t322 + t320 * t249) + t317 * (t320 * t141 - t249 * t322), t315 * (t188 * t322 + t320 * t204) + t317 * (t320 * t188 - t204 * t322), t315 * (t192 * t322 - t361) + t317 * (t320 * t192 + t363), t315 * (t189 * t322 + t320 * t343) + t317 * (t320 * t189 - t322 * t343), t315 * (t211 * t322 + t379) + t317 * (t211 * t320 - t258), t315 * (-pkin(7) * t124 + t113 * t322 - t320 * t95) + t317 * (-pkin(2) * t164 + pkin(7) * t125 + t320 * t113 + t322 * t95) - pkin(1) * t164 + qJ(2) * (-t124 * t315 + t125 * t317), t315 * (-pkin(7) * t126 - t320 * t101 + t116 * t322) + t317 * (-pkin(2) * t172 + pkin(7) * t127 + t101 * t322 + t320 * t116) - pkin(1) * t172 + qJ(2) * (-t126 * t315 + t127 * t317), t315 * (-pkin(7) * t114 + t322 * t59) + t317 * (pkin(7) * t115 + t320 * t59) + qJ(2) * (-t114 * t315 + t115 * t317) + (pkin(3) * t382 + t317 * t359 - pkin(1)) * t140, (t315 * (pkin(3) * t320 - pkin(8) * t322) + t317 * (-pkin(8) * t320 + t359) - pkin(1)) * t74 + t447 * (-t315 * (-t174 * t322 + t320 * t75) + t317 * (t320 * t174 + t322 * t75)), t330, -t513, t502, t450, t521, t449, t315 * (t30 * t322 - t320 * t31 - t506) + t317 * (t320 * t30 + t31 * t322 + t505) + t503, t315 * (t32 * t322 - t320 * t34 + t524) + t317 * (t320 * t32 + t322 * t34 - t523) - t525, t315 * (-t320 * t41 + t322 * t6 - t498) + t317 * (t320 * t6 + t322 * t41 + t494) + t488, t315 * (-pkin(7) * t12 + t3 * t322 - t320 * t5) + t317 * (-pkin(2) * t15 + pkin(7) * t13 + t320 * t3 + t322 * t5) - pkin(1) * t15 + qJ(2) * (-t12 * t315 + t13 * t317), t330, t502, t513, t449, -t521, t450, t315 * (t18 * t322 - t320 * t27 - t506) + t317 * (t320 * t18 + t27 * t322 + t505) + t503, t315 * (-t320 * t33 + t322 * t4 - t498) + t317 * (t320 * t4 + t322 * t33 + t494) + t488, t315 * (t17 * t322 - t320 * t25 - t524) + t317 * (t320 * t17 + t25 * t322 + t523) + t525, t315 * (-pkin(7) * t7 + t1 * t322 - t320 * t2) + t317 * (-pkin(2) * t9 + pkin(7) * t8 + t320 * t1 + t2 * t322) - pkin(1) * t9 + qJ(2) * (-t315 * t7 + t317 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t366, t367, -t300, -t293, 0, 0, 0, 0, 0, 0, t270, -t273, t250, -t266, 0, 0, 0, 0, 0, 0, t164, t172, t140, t74, 0, 0, 0, 0, 0, 0, t470, t81, t451, t15, 0, 0, 0, 0, 0, 0, t470, t451, -t81, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, t295 - t294, t296, -t275, -t347, qJDD(3), -t207, -t208, 0, 0, t233 * t319 + t282 * t383, t201 * t319 + t203 * t321, t255 * t321 + t463, t232 * t321 + t281 * t384, t254 * t319 + t395, (-t281 * t319 - t282 * t321) * t292, pkin(3) * t201 + pkin(8) * t165 - t400, pkin(3) * t205 + pkin(8) * t173 + t401, pkin(3) * t225 + pkin(8) * t142 + t75, -pkin(3) * t174 + pkin(8) * t75, t97, -t504, t490, t428, -t507, t431, t319 * t73 + t321 * t60 + t493, t319 * t80 + t321 * t66 - t516, t21 * t321 + t24 * t319 + t474, -pkin(3) * t110 + pkin(8) * t16 - qJ(5) * t413 + t26 * t321, t97, t490, t504, t431, t507, t428, t319 * t38 + t321 * t36 + t493, t19 * t321 + t20 * t319 + t474, t319 * t37 + t321 * t35 + t516, -pkin(3) * t61 + pkin(8) * t10 + t11 * t321 + t14 * t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t249, t204, -t251, t343, t262, -t118, -t119, 0, 0, t392, -t195, t439, -t392, t339, t262, t241 - t357 + t417, t240 + t352, t418, t421, t392, t439, t195, t262, -t339, -t392, -t194 * t247 + t241 - t329, t333, t436 + t444, t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, t441, t162, t110, 0, 0, 0, 0, 0, 0, t338, t162, -t441, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t440, t439, t443, t43;];
tauJ_reg  = t44;