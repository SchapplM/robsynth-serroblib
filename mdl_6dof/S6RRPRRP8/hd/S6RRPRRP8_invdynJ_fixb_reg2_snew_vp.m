% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 18:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:23:00
% EndTime: 2019-05-06 18:23:34
% DurationCPUTime: 17.68s
% Computational Cost: add. (69769->576), mult. (152793->764), div. (0->0), fcn. (112936->10), ass. (0->371)
t359 = cos(qJ(5));
t355 = sin(qJ(5));
t353 = sin(pkin(10));
t354 = cos(pkin(10));
t357 = sin(qJ(2));
t410 = qJD(1) * t357;
t322 = -t354 * qJD(2) + t353 * t410;
t323 = qJD(2) * t353 + t354 * t410;
t356 = sin(qJ(4));
t360 = cos(qJ(4));
t293 = t360 * t322 + t323 * t356;
t295 = -t322 * t356 + t323 * t360;
t264 = t359 * t293 + t295 * t355;
t266 = -t293 * t355 + t295 * t359;
t212 = t266 * t264;
t405 = qJD(1) * qJD(2);
t341 = t357 * t405;
t361 = cos(qJ(2));
t404 = t361 * qJDD(1);
t328 = -t341 + t404;
t324 = -qJDD(4) + t328;
t320 = -qJDD(5) + t324;
t472 = -t212 + t320;
t435 = t472 * t355;
t263 = t266 ^ 2;
t409 = qJD(1) * t361;
t339 = -qJD(4) + t409;
t335 = -qJD(5) + t339;
t455 = t335 ^ 2;
t471 = -t263 - t455;
t159 = t359 * t471 + t435;
t434 = t472 * t359;
t161 = -t355 * t471 + t434;
t108 = t159 * t360 + t161 * t356;
t99 = t159 * t356 - t161 * t360;
t72 = t108 * t354 - t353 * t99;
t538 = pkin(1) * t72;
t66 = t108 * t353 + t354 * t99;
t537 = qJ(3) * t66;
t536 = qJ(3) * t72;
t535 = t361 * t66;
t534 = -pkin(2) * t72 - pkin(3) * t108;
t458 = t264 ^ 2;
t248 = t458 - t455;
t174 = t248 * t355 - t434;
t178 = t248 * t359 + t435;
t117 = t174 * t360 + t178 * t356;
t120 = t174 * t356 - t178 * t360;
t345 = t357 * qJDD(1);
t398 = t361 * t405;
t327 = t345 + t398;
t302 = t353 * qJDD(2) + t354 * t327;
t392 = -t354 * qJDD(2) + t327 * t353;
t393 = t356 * t302 + t360 * t392;
t245 = -qJD(4) * t295 - t393;
t246 = -t293 * qJD(4) + t360 * t302 - t356 * t392;
t395 = -t359 * t245 + t355 * t246;
t407 = -qJD(5) - t335;
t370 = t266 * t407 - t395;
t533 = t357 * (t117 * t353 + t120 * t354) + t361 * t370;
t532 = pkin(8) * t99;
t530 = pkin(8) * t108;
t529 = t117 * t354 - t120 * t353;
t249 = -t263 + t455;
t199 = t320 + t212;
t433 = t199 * t355;
t510 = t359 * t249 - t433;
t195 = t359 * t199;
t511 = -t249 * t355 - t195;
t516 = -t356 * t510 + t360 * t511;
t517 = t356 * t511 + t360 * t510;
t527 = t353 * t516 + t354 * t517;
t526 = t357 * (-t353 * t517 + t354 * t516);
t468 = -t455 - t458;
t484 = t359 * t468 + t433;
t485 = t355 * t468 - t195;
t491 = t356 * t484 + t360 * t485;
t492 = -t356 * t485 + t360 * t484;
t507 = t353 * t492 + t354 * t491;
t525 = pkin(1) * t507;
t524 = pkin(4) * t159;
t523 = pkin(9) * t159;
t522 = pkin(9) * t161;
t506 = -t353 * t491 + t354 * t492;
t521 = qJ(3) * t506;
t520 = qJ(3) * t507;
t519 = t361 * t506;
t518 = -pkin(2) * t507 - pkin(3) * t491;
t514 = pkin(8) * t491;
t513 = pkin(8) * t492;
t379 = t355 * t245 + t359 * t246;
t182 = -qJD(5) * t264 + t379;
t251 = t264 * t335;
t144 = t182 - t251;
t463 = t144 * t355 + t359 * t370;
t464 = -t359 * t144 + t355 * t370;
t480 = t356 * t463 + t360 * t464;
t481 = -t356 * t464 + t360 * t463;
t494 = t353 * t481 + t354 * t480;
t512 = qJ(3) * t494;
t509 = -pkin(2) * t494 - pkin(3) * t480;
t184 = -t458 - t263;
t493 = -t353 * t480 + t354 * t481;
t508 = -pkin(2) * t184 + qJ(3) * t493;
t505 = pkin(7) * (t184 * t357 + t361 * t493) - pkin(1) * t494;
t502 = pkin(4) * t485;
t501 = pkin(8) * t480;
t500 = pkin(9) * t484;
t499 = pkin(9) * t485;
t363 = qJD(1) ^ 2;
t358 = sin(qJ(1));
t362 = cos(qJ(1));
t385 = g(1) * t362 + g(2) * t358;
t442 = qJDD(1) * pkin(7);
t317 = -pkin(1) * t363 - t385 + t442;
t382 = -pkin(2) * t361 - qJ(3) * t357;
t391 = t363 * t382 + t317;
t449 = t361 * g(3);
t454 = qJD(2) ^ 2;
t275 = -qJDD(2) * pkin(2) - t454 * qJ(3) + t357 * t391 + qJDD(3) + t449;
t303 = -pkin(3) * t409 - pkin(8) * t323;
t456 = t322 ^ 2;
t232 = t392 * pkin(3) - t456 * pkin(8) + t323 * t303 + t275;
t277 = -pkin(4) * t339 - pkin(9) * t295;
t457 = t293 ^ 2;
t155 = -t245 * pkin(4) - t457 * pkin(9) + t295 * t277 + t232;
t181 = -qJD(5) * t266 - t395;
t473 = t182 + t251;
t497 = -t181 * pkin(5) - qJ(6) * t473 + t155;
t496 = -pkin(4) * t184 + pkin(9) * t463;
t495 = -pkin(3) * t184 + pkin(8) * t481;
t452 = pkin(4) * t464;
t489 = pkin(9) * t464;
t470 = t263 - t458;
t486 = t361 * t470;
t417 = t335 * t359;
t403 = t264 * t417;
t375 = -t181 * t355 - t403;
t418 = t335 * t355;
t384 = t359 * t181 - t264 * t418;
t459 = t356 * t375 + t360 * t384;
t460 = -t356 * t384 + t360 * t375;
t483 = t353 * t460 + t354 * t459;
t372 = (t264 * t355 + t266 * t359) * t335;
t242 = t266 * t418;
t383 = -t242 + t403;
t461 = -t356 * t372 + t360 * t383;
t462 = t356 * t383 + t360 * t372;
t482 = t353 * t461 + t354 * t462;
t402 = t361 * t212;
t479 = t357 * (-t353 * t459 + t354 * t460) + t402;
t478 = t357 * (-t353 * t462 + t354 * t461) + t361 * t320;
t420 = t322 * t323;
t373 = -t328 - t420;
t477 = t353 * t373;
t476 = t354 * t373;
t421 = t295 * t293;
t371 = -t324 - t421;
t475 = t356 * t371;
t474 = t360 * t371;
t282 = t293 * t339;
t229 = -t282 + t246;
t469 = t282 + t246;
t399 = t322 * t409;
t284 = -t302 + t399;
t307 = t323 * t409;
t286 = -t392 - t307;
t292 = t295 ^ 2;
t321 = t323 ^ 2;
t337 = t339 ^ 2;
t396 = t358 * g(1) - t362 * g(2);
t316 = qJDD(1) * pkin(1) + t363 * pkin(7) + t396;
t381 = t327 + t398;
t272 = -t381 * qJ(3) + (-t328 + t341) * pkin(2) - t316;
t450 = t357 * g(3);
t276 = -pkin(2) * t454 + qJDD(2) * qJ(3) + t361 * t391 - t450;
t230 = 0.2e1 * qJD(3) * t323 - t354 * t272 + t353 * t276;
t188 = t373 * pkin(3) + pkin(8) * t284 - t230;
t231 = -0.2e1 * qJD(3) * t322 + t353 * t272 + t354 * t276;
t190 = -pkin(3) * t456 - pkin(8) * t392 + t303 * t409 + t231;
t124 = -t360 * t188 + t356 * t190;
t104 = t371 * pkin(4) - pkin(9) * t229 - t124;
t125 = t356 * t188 + t360 * t190;
t107 = -pkin(4) * t457 + pkin(9) * t245 + t277 * t339 + t125;
t69 = -t359 * t104 + t107 * t355;
t70 = t355 * t104 + t359 * t107;
t37 = t355 * t70 - t359 * t69;
t453 = pkin(4) * t37;
t451 = pkin(5) * t359;
t408 = qJD(6) * t335;
t330 = -0.2e1 * t408;
t208 = pkin(5) * t264 - qJ(6) * t266;
t389 = -t320 * qJ(6) - t264 * t208 + t70;
t374 = -pkin(5) * t455 + t389;
t49 = t330 + t374;
t51 = t320 * pkin(5) - qJ(6) * t455 + t208 * t266 + qJDD(6) + t69;
t448 = -pkin(5) * t51 + qJ(6) * t49;
t81 = -t124 * t360 + t125 * t356;
t447 = t353 * t81;
t446 = t354 * t81;
t445 = t356 * t37;
t444 = t360 * t37;
t443 = qJ(6) * t359;
t406 = qJD(5) - t335;
t140 = t266 * t406 + t395;
t441 = t140 * t355;
t440 = t140 * t359;
t147 = t264 * t406 - t379;
t438 = t147 * t355;
t437 = t155 * t355;
t436 = t155 * t359;
t432 = t232 * t356;
t431 = t232 * t360;
t254 = t324 - t421;
t429 = t254 * t356;
t428 = t254 * t360;
t427 = t266 * t335;
t426 = t275 * t353;
t425 = t275 * t354;
t287 = t328 - t420;
t424 = t287 * t353;
t423 = t287 * t354;
t416 = t339 * t356;
t415 = t339 * t360;
t338 = t361 * t363 * t357;
t414 = t357 * (qJDD(2) + t338);
t412 = t361 * (-t338 + qJDD(2));
t411 = -pkin(5) * t144 + qJ(6) * t370;
t401 = t361 * t421;
t400 = t361 * t420;
t397 = -qJ(6) * t355 - pkin(4);
t38 = t355 * t69 + t359 * t70;
t82 = t124 * t356 + t360 * t125;
t169 = t230 * t353 + t354 * t231;
t299 = t357 * t317 + t449;
t300 = t317 * t361 - t450;
t394 = t357 * t299 + t361 * t300;
t27 = t355 * t49 - t359 * t51;
t390 = pkin(4) * t27 + t448;
t388 = t411 + t452;
t387 = -t70 + t524;
t135 = t182 * t355 - t266 * t417;
t136 = t182 * t359 + t242;
t85 = t135 * t360 + t136 * t356;
t88 = -t135 * t356 + t136 * t360;
t386 = t357 * (-t353 * t85 + t354 * t88) - t402;
t380 = t230 * t354 - t231 * t353;
t378 = -pkin(1) + t382;
t376 = -t69 + t502;
t369 = (-qJD(4) - t339) * t295 - t393;
t368 = -pkin(5) * t471 - qJ(6) * t472 + t374;
t367 = t368 - t524;
t366 = -pkin(5) * t199 + qJ(6) * t468 - t51;
t365 = t366 + t502;
t364 = 0.2e1 * qJD(6) * t266 - t497;
t351 = t361 ^ 2;
t350 = t357 ^ 2;
t347 = t351 * t363;
t346 = t350 * t363;
t329 = -0.2e1 * t341 + t404;
t326 = t345 + 0.2e1 * t398;
t318 = t361 * t328;
t306 = -t321 - t347;
t305 = -t321 + t347;
t304 = -t347 + t456;
t296 = -t347 - t456;
t285 = -t307 + t392;
t283 = -t399 - t302;
t280 = -t321 - t456;
t279 = -t292 + t337;
t278 = -t337 + t457;
t274 = -t292 - t337;
t269 = -t306 * t353 + t423;
t268 = t306 * t354 + t424;
t267 = t292 - t457;
t260 = -t337 - t457;
t258 = t296 * t354 - t477;
t257 = t296 * t353 + t476;
t244 = -t284 * t353 + t286 * t354;
t239 = (t293 * t360 - t295 * t356) * t339;
t238 = (t293 * t356 + t295 * t360) * t339;
t233 = -t292 - t457;
t224 = (qJD(4) - t339) * t295 + t393;
t222 = t278 * t360 + t429;
t221 = -t279 * t356 + t474;
t220 = t278 * t356 - t428;
t219 = t279 * t360 + t475;
t218 = t246 * t360 + t295 * t416;
t217 = t246 * t356 - t295 * t415;
t216 = -t245 * t356 - t293 * t415;
t215 = t245 * t360 - t293 * t416;
t214 = -t274 * t356 + t428;
t213 = t274 * t360 + t429;
t204 = t260 * t360 - t475;
t203 = t260 * t356 + t474;
t167 = t229 * t356 + t360 * t369;
t166 = -t224 * t360 - t356 * t469;
t165 = -t229 * t360 + t356 * t369;
t164 = -t224 * t356 + t360 * t469;
t163 = -pkin(8) * t213 + t431;
t158 = -t213 * t353 + t214 * t354;
t157 = t213 * t354 + t214 * t353;
t156 = -pkin(8) * t203 + t432;
t150 = -t203 * t353 + t204 * t354;
t149 = t203 * t354 + t204 * t353;
t145 = t264 * t407 + t379;
t141 = -t181 - t427;
t126 = -pkin(3) * t469 + pkin(8) * t214 + t432;
t122 = -pkin(3) * t224 + pkin(8) * t204 - t431;
t113 = -t165 * t353 + t167 * t354;
t112 = t165 * t354 + t167 * t353;
t105 = t436 - t523;
t101 = t437 - t499;
t94 = -t355 * t473 - t440;
t93 = -t438 + t440;
t90 = t359 * t473 - t441;
t89 = t147 * t359 + t441;
t80 = -pkin(4) * t473 + t437 + t522;
t79 = -pkin(3) * t232 + pkin(8) * t82;
t78 = -pkin(4) * t140 - t436 + t500;
t77 = -pkin(8) * t165 - t81;
t76 = (-pkin(5) * t335 - 0.2e1 * qJD(6)) * t266 + t497;
t71 = -pkin(3) * t233 + pkin(8) * t167 + t82;
t61 = -t356 * t90 + t360 * t94;
t60 = -t356 * t89 + t360 * t93;
t57 = t356 * t94 + t360 * t90;
t56 = t356 * t93 + t360 * t89;
t55 = t353 * t88 + t354 * t85;
t53 = (-t141 + t427) * pkin(5) + t364;
t52 = pkin(5) * t427 - qJ(6) * t147 + t364;
t47 = t354 * t82 - t447;
t46 = t353 * t82 + t446;
t45 = -qJ(6) * t184 + t51;
t44 = t330 + (-t184 - t455) * pkin(5) + t389;
t43 = -t141 * t443 - t355 * t53 - t499;
t42 = t105 * t360 - t356 * t80 - t530;
t41 = pkin(5) * t438 + t359 * t52 + t523;
t40 = t101 * t360 - t356 * t78 - t514;
t39 = t141 * t397 + t359 * t53 + t500;
t36 = -pkin(3) * t473 + t105 * t356 + t360 * t80 - t532;
t35 = -t522 + t355 * t52 + (-pkin(4) - t451) * t147;
t34 = -pkin(3) * t140 + t101 * t356 + t360 * t78 + t513;
t33 = -pkin(4) * t155 + pkin(9) * t38;
t28 = t355 * t51 + t359 * t49;
t26 = -t37 - t489;
t25 = t38 + t496;
t24 = -t355 * t44 + t359 * t45 - t489;
t23 = t355 * t45 + t359 * t44 + t496;
t22 = t360 * t38 - t445;
t21 = t356 * t38 + t444;
t20 = -t356 * t39 + t360 * t43 - t514;
t19 = -pkin(9) * t27 + (pkin(5) * t355 - t443) * t76;
t18 = -t35 * t356 + t360 * t41 + t530;
t17 = -pkin(3) * t141 + t356 * t43 + t360 * t39 + t513;
t16 = -pkin(3) * t147 + t35 * t360 + t356 * t41 + t532;
t15 = -t27 * t356 + t28 * t360;
t14 = t27 * t360 + t28 * t356;
t13 = pkin(9) * t28 + (t397 - t451) * t76;
t12 = -t25 * t356 + t26 * t360 - t501;
t11 = t25 * t360 + t26 * t356 + t495;
t10 = -t21 * t353 + t22 * t354;
t9 = t21 * t354 + t22 * t353;
t8 = -t23 * t356 + t24 * t360 - t501;
t7 = -pkin(8) * t21 - pkin(9) * t444 - t33 * t356;
t6 = t23 * t360 + t24 * t356 + t495;
t5 = -pkin(3) * t155 + pkin(8) * t22 - pkin(9) * t445 + t33 * t360;
t4 = -t14 * t353 + t15 * t354;
t3 = t14 * t354 + t15 * t353;
t2 = -pkin(8) * t14 - t13 * t356 + t19 * t360;
t1 = -pkin(3) * t76 + pkin(8) * t15 + t13 * t360 + t19 * t356;
t29 = [0, 0, 0, 0, 0, qJDD(1), t396, t385, 0, 0, t381 * t357, t326 * t361 + t329 * t357, t414 + t361 * (-t346 + t454), -t341 * t361 + t318, t357 * (t347 - t454) + t412, 0, t361 * t316 + pkin(1) * t329 + pkin(7) * (t361 * (-t347 - t454) - t414), -t357 * t316 - pkin(1) * t326 + pkin(7) * (-t412 - t357 * (-t346 - t454)), pkin(1) * (t346 + t347) + (t350 + t351) * t442 + t394, pkin(1) * t316 + pkin(7) * t394, t357 * (t302 * t354 + t307 * t353) - t400, t357 * (t283 * t353 - t285 * t354) + t361 * (-t321 + t456), t357 * (-t305 * t353 + t476) + t361 * t284, t357 * (t353 * t392 - t354 * t399) + t400, t357 * (t304 * t354 + t424) - t361 * t286, t318 + t357 * (t322 * t354 - t323 * t353) * t409, t357 * (-qJ(3) * t257 + t426) + t361 * (-pkin(2) * t257 + t230) - pkin(1) * t257 + pkin(7) * (t258 * t361 + t285 * t357), t357 * (-qJ(3) * t268 + t425) + t361 * (-pkin(2) * t268 + t231) - pkin(1) * t268 + pkin(7) * (t269 * t361 - t283 * t357), t357 * t380 + pkin(7) * (t244 * t361 + t280 * t357) + t378 * (t284 * t354 + t286 * t353), pkin(7) * (t169 * t361 + t275 * t357) - t378 * t380, t357 * (-t217 * t353 + t218 * t354) - t401, t357 * (-t164 * t353 + t166 * t354) - t361 * t267, t357 * (-t219 * t353 + t221 * t354) - t361 * t229, t357 * (-t215 * t353 + t216 * t354) + t401, t357 * (-t220 * t353 + t222 * t354) - t361 * t369, t357 * (-t238 * t353 + t239 * t354) + t361 * t324, t357 * (-qJ(3) * t149 - t122 * t353 + t156 * t354) + t361 * (-pkin(2) * t149 - pkin(3) * t203 + t124) - pkin(1) * t149 + pkin(7) * (t150 * t361 + t224 * t357), t357 * (-qJ(3) * t157 - t126 * t353 + t163 * t354) + t361 * (-pkin(2) * t157 - pkin(3) * t213 + t125) - pkin(1) * t157 + pkin(7) * (t158 * t361 + t357 * t469), t357 * (-qJ(3) * t112 - t353 * t71 + t354 * t77) + t361 * (-pkin(2) * t112 - pkin(3) * t165) - pkin(1) * t112 + pkin(7) * (t113 * t361 + t233 * t357), t357 * (-pkin(8) * t446 - qJ(3) * t46 - t353 * t79) + t361 * (-pkin(2) * t46 - pkin(3) * t81) - pkin(1) * t46 + pkin(7) * (t232 * t357 + t361 * t47), t386, t357 * (-t353 * t57 + t354 * t61) - t486, -t361 * t144 + t526, t479, -t533, t478, t357 * (-t34 * t353 + t354 * t40 - t520) + t361 * (-t376 + t518) - t525 + pkin(7) * (t140 * t357 + t519), t357 * (-t353 * t36 + t354 * t42 - t536) + t361 * (-t387 + t534) - t538 + pkin(7) * (t357 * t473 - t535), t357 * (-t11 * t353 + t12 * t354 - t512) + t361 * (-t452 + t509) + t505, t357 * (-qJ(3) * t9 - t353 * t5 + t354 * t7) + t361 * (-pkin(2) * t9 - pkin(3) * t21 - t453) - pkin(1) * t9 + pkin(7) * (t10 * t361 + t155 * t357), t386, -t361 * t145 + t526, t357 * (-t353 * t56 + t354 * t60) + t486, t478, t533, t479, t357 * (-t17 * t353 + t20 * t354 - t520) + t361 * (-t365 + t518) - t525 + pkin(7) * (t141 * t357 + t519), t357 * (-t353 * t6 + t354 * t8 - t512) + t361 * (-t388 + t509) + t505, t357 * (-t16 * t353 + t18 * t354 + t536) + t361 * (-t367 + 0.2e1 * t408 - t534) + t538 + pkin(7) * (t147 * t357 + t535), t357 * (-qJ(3) * t3 - t1 * t353 + t2 * t354) + t361 * (-pkin(2) * t3 - pkin(3) * t14 - t390) - pkin(1) * t3 + pkin(7) * (t357 * t76 + t361 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t338, t346 - t347, t345, t338, t404, qJDD(2), -t299, -t300, 0, 0, t302 * t353 - t307 * t354, -t283 * t354 - t285 * t353, t305 * t354 + t477, -t353 * t399 - t354 * t392, t304 * t353 - t423, (t322 * t353 + t323 * t354) * t409, -pkin(2) * t285 + qJ(3) * t258 - t425, pkin(2) * t283 + qJ(3) * t269 + t426, -pkin(2) * t280 + qJ(3) * t244 + t169, -pkin(2) * t275 + qJ(3) * t169, t217 * t354 + t218 * t353, t164 * t354 + t166 * t353, t219 * t354 + t221 * t353, t215 * t354 + t216 * t353, t220 * t354 + t222 * t353, t238 * t354 + t239 * t353, -pkin(2) * t224 + qJ(3) * t150 + t122 * t354 + t156 * t353, -pkin(2) * t469 + qJ(3) * t158 + t126 * t354 + t163 * t353, -pkin(2) * t233 + qJ(3) * t113 + t353 * t77 + t354 * t71, -pkin(2) * t232 - pkin(8) * t447 + qJ(3) * t47 + t354 * t79, t55, t353 * t61 + t354 * t57, t527, t483, t529, t482, -pkin(2) * t140 + t34 * t354 + t353 * t40 + t521, -pkin(2) * t473 + t353 * t42 + t354 * t36 - t537, t11 * t354 + t12 * t353 + t508, -pkin(2) * t155 + qJ(3) * t10 + t353 * t7 + t354 * t5, t55, t527, t353 * t60 + t354 * t56, t482, -t529, t483, -pkin(2) * t141 + t17 * t354 + t20 * t353 + t521, t353 * t8 + t354 * t6 + t508, -pkin(2) * t147 + t16 * t354 + t18 * t353 + t537, -pkin(2) * t76 + qJ(3) * t4 + t1 * t354 + t2 * t353; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, -t283, t280, t275, 0, 0, 0, 0, 0, 0, t224, t469, t233, t232, 0, 0, 0, 0, 0, 0, t140, t473, t184, t155, 0, 0, 0, 0, 0, 0, t141, t184, t147, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t421, t267, t229, -t421, t369, -t324, -t124, -t125, 0, 0, t212, t470, t144, -t212, t370, -t320, t376, t387, t452, t453, t212, t145, -t470, -t320, -t370, -t212, t365, t388, t330 + t367, t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t470, t144, -t212, t370, -t320, -t69, -t70, 0, 0, t212, t145, -t470, -t320, -t370, -t212, t366, t411, t330 + t368, t448; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, t144, t471, t51;];
tauJ_reg  = t29;
