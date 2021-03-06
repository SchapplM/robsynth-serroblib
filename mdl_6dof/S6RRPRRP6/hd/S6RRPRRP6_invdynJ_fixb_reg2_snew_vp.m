% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPRRP6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 18:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPRRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:05:55
% EndTime: 2019-05-06 18:06:56
% DurationCPUTime: 25.04s
% Computational Cost: add. (78346->671), mult. (205911->940), div. (0->0), fcn. (164440->12), ass. (0->406)
t366 = sin(qJ(2));
t369 = cos(qJ(2));
t364 = sin(qJ(5));
t360 = sin(pkin(11));
t362 = cos(pkin(11));
t361 = sin(pkin(6));
t440 = qJD(1) * t369;
t417 = t361 * t440;
t441 = qJD(1) * t366;
t418 = t361 * t441;
t332 = t360 * t417 + t362 * t418;
t363 = cos(pkin(6));
t442 = qJD(1) * t363;
t355 = qJD(2) + t442;
t365 = sin(qJ(4));
t368 = cos(qJ(4));
t313 = t332 * t368 + t355 * t365;
t428 = qJDD(1) * t369;
t431 = qJD(1) * qJD(2);
t382 = (-t366 * t431 + t428) * t361;
t429 = qJDD(1) * t366;
t383 = (t369 * t431 + t429) * t361;
t303 = t360 * t382 + t362 * t383;
t354 = t363 * qJDD(1) + qJDD(2);
t410 = t303 * t365 - t368 * t354;
t267 = -qJD(4) * t313 - t410;
t266 = qJDD(5) - t267;
t330 = t360 * t418 - t362 * t417;
t326 = qJD(4) + t330;
t367 = cos(qJ(5));
t286 = t313 * t364 - t367 * t326;
t288 = t313 * t367 + t326 * t364;
t461 = t288 * t286;
t198 = -t461 - t266;
t472 = t198 * t367;
t285 = t288 ^ 2;
t311 = t332 * t365 - t368 * t355;
t308 = qJD(5) + t311;
t502 = t308 ^ 2;
t524 = -t285 - t502;
t146 = -t364 * t524 + t472;
t393 = -t303 * t368 - t354 * t365;
t268 = -qJD(4) * t311 - t393;
t380 = t360 * t383 - t362 * t382;
t378 = -qJDD(4) - t380;
t372 = -t286 * qJD(5) + t367 * t268 - t364 * t378;
t462 = t286 * t308;
t517 = -t462 + t372;
t100 = t146 * t368 + t365 * t517;
t473 = t198 * t364;
t144 = t367 * t524 + t473;
t76 = t100 * t360 - t144 * t362;
t78 = t100 * t362 + t144 * t360;
t610 = pkin(1) * (t366 * t78 + t369 * t76);
t98 = t146 * t365 - t368 * t517;
t609 = pkin(1) * t98 + pkin(8) * (t366 * t76 - t369 * t78);
t608 = qJ(3) * t76;
t607 = -pkin(2) * t76 + pkin(3) * t144 - pkin(9) * t100;
t606 = pkin(2) * t98 - qJ(3) * t78;
t503 = t286 ^ 2;
t253 = t503 - t502;
t158 = -t253 * t367 - t473;
t257 = t308 * t288;
t411 = -t268 * t364 - t367 * t378;
t390 = qJD(5) * t288 - t411;
t173 = -t257 + t390;
t104 = t158 * t365 - t173 * t368;
t108 = t158 * t368 + t173 * t365;
t154 = -t253 * t364 + t472;
t605 = t363 * t104 + (t366 * (t108 * t362 + t154 * t360) + t369 * (t108 * t360 - t154 * t362)) * t361;
t526 = t257 + t390;
t115 = -t526 * t364 + t367 * t517;
t474 = t517 * t364;
t117 = t526 * t367 + t474;
t523 = t285 - t503;
t86 = t117 * t365 + t368 * t523;
t88 = t117 * t368 - t365 * t523;
t602 = t361 * (t366 * (-t115 * t360 + t362 * t88) + t369 * (t115 * t362 + t360 * t88)) + t363 * t86;
t598 = pkin(3) * t98;
t597 = pkin(9) * t98;
t587 = pkin(4) * t144;
t586 = pkin(10) * t144;
t585 = pkin(10) * t146;
t518 = -t461 + t266;
t470 = t518 * t367;
t515 = -t502 - t503;
t531 = t364 * t515 + t470;
t471 = t518 * t364;
t530 = t367 * t515 - t471;
t548 = t365 * t526 + t368 * t530;
t565 = t360 * t548 - t362 * t531;
t580 = qJ(3) * t565;
t579 = pkin(2) * t565 - pkin(3) * t531 + pkin(9) * t548;
t549 = t365 * t530 - t368 * t526;
t564 = t360 * t531 + t362 * t548;
t578 = -pkin(2) * t549 + qJ(3) * t564;
t577 = pkin(1) * (t366 * t564 + t369 * t565);
t576 = pkin(8) * (-t366 * t565 + t369 * t564) - pkin(1) * t549;
t254 = -t285 + t502;
t551 = t254 * t367 + t471;
t516 = t462 + t372;
t550 = -t254 * t364 + t470;
t562 = t365 * t516 + t368 * t550;
t563 = t365 * t550 - t368 * t516;
t575 = t363 * t563 + (t366 * (t360 * t551 + t362 * t562) + t369 * (t360 * t562 - t362 * t551)) * t361;
t572 = pkin(3) * t549;
t571 = pkin(9) * t549;
t560 = pkin(4) * t531;
t559 = pkin(10) * t530;
t558 = pkin(10) * t531;
t434 = -qJD(2) - t355;
t522 = t285 + t503;
t547 = pkin(4) * t522;
t546 = qJ(6) * t517;
t302 = t332 * t330;
t521 = -t302 + t354;
t545 = t360 * t521;
t544 = t362 * t521;
t541 = t365 * t522;
t274 = t313 * t311;
t525 = -t274 - t378;
t539 = t365 * t525;
t535 = t368 * t522;
t533 = t368 * t525;
t458 = t308 * t364;
t167 = t286 * t458 - t367 * t390;
t457 = t308 * t367;
t425 = t286 * t457;
t391 = t364 * t390 + t425;
t427 = t365 * t461;
t506 = t368 * t391 - t427;
t426 = t368 * t461;
t507 = t365 * t391 + t426;
t529 = t363 * t507 + (t366 * (t167 * t360 + t362 * t506) + t369 * (-t167 * t362 + t360 * t506)) * t361;
t250 = t288 * t458;
t400 = t250 - t425;
t505 = t266 * t365 + t368 * t400;
t509 = -t368 * t266 + t365 * t400;
t527 = (t286 * t364 + t288 * t367) * t308;
t528 = t363 * t509 + (t366 * (-t360 * t527 + t362 * t505) + t369 * (t360 * t505 + t362 * t527)) * t361;
t291 = t326 * t311;
t226 = t268 - t291;
t319 = t355 * t330;
t520 = -t319 + t303;
t358 = t361 ^ 2;
t498 = qJD(1) ^ 2;
t436 = t358 * t498;
t414 = t369 * t436;
t350 = t366 * t414;
t337 = t354 + t350;
t454 = t355 * t332;
t275 = t380 + t454;
t494 = sin(qJ(1));
t495 = cos(qJ(1));
t388 = g(1) * t495 + g(2) * t494;
t430 = qJDD(1) * t361;
t336 = -pkin(1) * t498 + pkin(8) * t430 - t388;
t451 = t361 * t366;
t356 = g(3) * t451;
t387 = g(1) * t494 - g(2) * t495;
t435 = t361 * t498;
t376 = qJDD(1) * pkin(1) + pkin(8) * t435 + t387;
t373 = t363 * t376;
t407 = qJD(1) * (-qJD(2) + t355);
t500 = t355 ^ 2;
t251 = -pkin(2) * t500 - t356 + (qJ(3) * t361 * t407 + t373) * t366 + (-pkin(2) * t414 + qJ(3) * t430 + t336) * t369;
t409 = -t366 * t336 + t369 * t373;
t371 = t337 * pkin(2) + (-t369 * g(3) + (t369 * t407 - t429) * qJ(3)) * t361 + t409;
t191 = -0.2e1 * qJD(3) * t330 + t362 * t251 + t360 * t371;
t450 = t361 * t369;
t292 = g(3) * t450 - t409;
t293 = t369 * t336 + t366 * t373 - t356;
t519 = t366 * t292 + t369 * t293;
t499 = t369 ^ 2;
t416 = t499 * t436;
t514 = qJ(3) * t416 - qJDD(3);
t233 = pkin(5) * t286 - qJ(6) * t288;
t299 = pkin(3) * t330 - pkin(9) * t332;
t412 = t251 * t360 - t362 * t371;
t497 = 0.2e1 * qJD(3);
t161 = -t354 * pkin(3) - t500 * pkin(9) + (t497 + t299) * t332 + t412;
t112 = -t226 * pkin(10) + (t313 * t326 - t267) * pkin(4) + t161;
t162 = -pkin(3) * t500 + pkin(9) * t354 - t299 * t330 + t191;
t487 = t363 * g(3);
t320 = t361 * t376 + t487;
t370 = (pkin(2) * t355 - qJ(3) * t418) * t418 - pkin(2) * t382 - t320 - t520 * pkin(9) + t275 * pkin(3) - t514;
t123 = t368 * t162 + t365 * t370;
t271 = pkin(4) * t311 - pkin(10) * t313;
t501 = t326 ^ 2;
t93 = -pkin(4) * t501 - pkin(10) * t378 - t311 * t271 + t123;
t65 = t364 * t112 + t367 * t93;
t406 = -t266 * qJ(6) + t286 * t233 - t65;
t513 = -pkin(5) * (t524 + t502) - qJ(6) * t198 - t406;
t437 = qJD(6) * t308;
t304 = 0.2e1 * t437;
t397 = t304 - t406;
t51 = -pkin(5) * t502 + t397;
t64 = -t367 * t112 + t364 * t93;
t54 = -t266 * pkin(5) - qJ(6) * t502 + t233 * t288 + qJDD(6) + t64;
t34 = t364 * t54 + t367 * t51;
t413 = qJ(6) * t364 + pkin(4);
t489 = pkin(5) * t367;
t122 = t162 * t365 - t368 * t370;
t92 = t378 * pkin(4) - t501 * pkin(10) + t271 * t313 + t122;
t379 = t390 * pkin(5) - t546 + t92;
t62 = (pkin(5) * t308 - 0.2e1 * qJD(6)) * t288 + t379;
t512 = -(t413 + t489) * t62 + pkin(10) * t34;
t377 = 0.2e1 * qJD(6) * t288 - t379;
t53 = (-t526 - t257) * pkin(5) + t377;
t511 = t367 * t53 - t413 * t526 + t559;
t52 = -pkin(5) * t257 + t377 + t546;
t510 = -t585 + t517 * (pkin(4) + t489) + t364 * t52;
t223 = (qJD(4) - t326) * t313 + t410;
t508 = (t442 + t434) * t358;
t170 = t288 * t457 + t364 * t372;
t171 = t367 * t372 - t250;
t401 = t368 * t171 + t427;
t402 = t365 * t171 - t426;
t504 = t363 * t402 + (t366 * (t170 * t360 + t362 * t401) + t369 * (-t362 * t170 + t360 * t401)) * t361;
t309 = t311 ^ 2;
t310 = t313 ^ 2;
t328 = t330 ^ 2;
t329 = t332 ^ 2;
t359 = t366 ^ 2;
t190 = t332 * t497 + t412;
t126 = -t190 * t362 + t191 * t360;
t493 = pkin(2) * t126;
t280 = t319 + t303;
t375 = -t380 + t454;
t229 = -t280 * t362 + t360 * t375;
t492 = pkin(2) * t229;
t491 = pkin(3) * t360;
t490 = pkin(4) * t365;
t488 = pkin(8) * t361;
t486 = t364 * t92;
t485 = t367 * t92;
t484 = qJ(6) * t367;
t483 = t161 * t365;
t482 = t161 * t368;
t476 = t516 * t364;
t475 = t516 * t367;
t242 = t274 - t378;
t469 = t242 * t365;
t468 = t242 * t368;
t270 = t487 + ((t369 * pkin(2) + pkin(1)) * qJDD(1) + ((qJ(3) * t359 + pkin(8)) * t361 * qJD(1) + t434 * t366 * pkin(2)) * qJD(1) + t387) * t361 + t514;
t464 = t270 * t360;
t463 = t270 * t362;
t296 = t302 + t354;
t460 = t296 * t360;
t459 = t296 * t362;
t456 = t326 * t365;
t455 = t326 * t368;
t453 = t355 * t360;
t452 = t355 * t362;
t449 = t366 * t126;
t447 = t366 * t337;
t338 = -t350 + t354;
t445 = t369 * t338;
t432 = qJD(4) + t326;
t424 = t360 * t274;
t423 = t362 * t274;
t422 = t363 * t302;
t421 = -pkin(3) * t362 - pkin(2);
t420 = -pkin(4) * t368 - pkin(3);
t415 = t359 * t436;
t39 = t364 * t64 + t367 * t65;
t73 = t122 * t365 + t368 * t123;
t127 = t190 * t360 + t362 * t191;
t404 = -pkin(4) * t92 + pkin(10) * t39;
t403 = -pkin(5) * t54 + qJ(6) * t51;
t396 = -pkin(5) * t516 - qJ(6) * t173;
t38 = t364 * t65 - t367 * t64;
t72 = -t122 * t368 + t123 * t365;
t386 = -pkin(4) * t526 - t485 + t559;
t385 = -pkin(4) * t517 + t486 + t585;
t175 = (-qJD(5) + t308) * t288 + t411;
t120 = t175 * t367 + t476;
t384 = pkin(10) * t120 + t39 + t547;
t118 = -t173 * t367 + t476;
t47 = (t522 - t502) * pkin(5) + t397;
t48 = qJ(6) * t522 + t54;
t381 = pkin(10) * t118 + t364 * t48 + t367 * t47 + t547;
t374 = pkin(5) * t518 + qJ(6) * t515 - t54;
t344 = t363 * t354;
t343 = t355 * t417;
t342 = t355 * t418;
t340 = (t359 - t499) * t436;
t339 = -t500 - t416;
t327 = -t500 - t415;
t318 = -t329 + t500;
t317 = t328 - t500;
t316 = -t342 + t382;
t315 = t342 + t382;
t314 = -t343 + t383;
t307 = -t329 - t500;
t301 = t329 - t328;
t294 = -t500 - t328;
t290 = -t310 + t501;
t289 = t309 - t501;
t282 = -t328 - t329;
t273 = t310 - t309;
t269 = -t310 - t501;
t265 = -t501 - t309;
t264 = -t307 * t360 - t459;
t263 = t307 * t362 - t460;
t249 = t309 + t310;
t245 = t294 * t362 - t545;
t244 = t294 * t360 + t544;
t240 = (-t311 * t368 + t313 * t365) * t326;
t239 = (-t311 * t365 - t313 * t368) * t326;
t230 = t280 * t360 + t362 * t375;
t228 = t311 * t432 + t393;
t227 = t268 + t291;
t224 = -t313 * t432 - t410;
t222 = t268 * t368 - t313 * t456;
t221 = t268 * t365 + t313 * t455;
t220 = -t267 * t365 + t311 * t455;
t219 = t267 * t368 + t311 * t456;
t216 = t289 * t368 - t469;
t215 = -t290 * t365 + t533;
t214 = t289 * t365 + t468;
t213 = t290 * t368 + t539;
t209 = -t269 * t365 - t468;
t208 = t269 * t368 - t469;
t200 = t265 * t368 - t539;
t199 = t265 * t365 + t533;
t188 = -t223 * t368 + t227 * t365;
t187 = t224 * t368 - t226 * t365;
t186 = -t223 * t365 - t227 * t368;
t185 = t224 * t365 + t226 * t368;
t151 = pkin(2) * t263 - t191;
t150 = pkin(2) * t244 - t190;
t149 = t209 * t362 - t228 * t360;
t148 = t209 * t360 + t228 * t362;
t143 = t200 * t362 - t224 * t360;
t142 = t200 * t360 + t224 * t362;
t137 = t188 * t362 - t249 * t360;
t136 = t188 * t360 + t249 * t362;
t125 = -pkin(9) * t208 + t482;
t124 = -pkin(9) * t199 + t483;
t116 = t175 * t364 - t475;
t114 = -t173 * t364 - t475;
t90 = -pkin(3) * t208 + t123;
t85 = -pkin(3) * t199 + t122;
t84 = t120 * t368 - t541;
t83 = t118 * t368 - t541;
t82 = t120 * t365 + t535;
t81 = t118 * t365 + t535;
t80 = pkin(2) * t148 + pkin(3) * t228 + pkin(9) * t209 + t483;
t75 = pkin(2) * t142 + pkin(3) * t224 + pkin(9) * t200 - t482;
t74 = t485 - t586;
t67 = t486 - t558;
t66 = -pkin(4) * t114 - t396;
t61 = -pkin(9) * t186 - t72;
t60 = t161 * t360 + t362 * t73;
t59 = -t161 * t362 + t360 * t73;
t58 = t116 * t360 + t362 * t84;
t57 = t114 * t360 + t362 * t83;
t56 = -t116 * t362 + t360 * t84;
t55 = -t114 * t362 + t360 * t83;
t50 = t65 - t587;
t49 = t64 - t560;
t46 = pkin(2) * t136 + pkin(3) * t249 + pkin(9) * t188 + t73;
t45 = -t385 - t598;
t44 = -t386 - t572;
t43 = -t364 * t53 - t484 * t526 - t558;
t42 = -t374 - t560;
t41 = -pkin(5) * t474 + t367 * t52 + t586;
t40 = -0.2e1 * t437 - t513 + t587;
t37 = pkin(2) * t59 - pkin(3) * t161 + pkin(9) * t73;
t36 = -pkin(10) * t116 - t38;
t35 = -t511 - t572;
t33 = t364 * t51 - t367 * t54;
t32 = -t510 + t598;
t31 = t365 * t92 + t368 * t39;
t30 = t365 * t39 - t368 * t92;
t29 = -t365 * t50 + t368 * t74 - t597;
t28 = -t365 * t49 + t368 * t67 - t571;
t27 = -pkin(3) * t82 - t384;
t26 = -pkin(10) * t114 - t364 * t47 + t367 * t48;
t25 = -pkin(9) * t82 + t116 * t490 + t36 * t368;
t24 = t34 * t368 + t365 * t62;
t23 = t34 * t365 - t368 * t62;
t22 = t365 * t74 + t368 * t50 - t607;
t21 = t365 * t67 + t368 * t49 + t579;
t20 = -pkin(3) * t81 - t381;
t19 = -t365 * t42 + t368 * t43 - t571;
t18 = -t365 * t40 + t368 * t41 + t597;
t17 = t31 * t362 + t360 * t38;
t16 = t31 * t360 - t362 * t38;
t15 = -pkin(10) * t33 + (pkin(5) * t364 - t484) * t62;
t14 = -pkin(9) * t81 + t26 * t368 - t365 * t66;
t13 = -pkin(4) * t33 - t403;
t12 = pkin(2) * t56 + pkin(9) * t84 + t116 * t420 + t36 * t365;
t11 = -pkin(3) * t30 - t404;
t10 = t365 * t43 + t368 * t42 + t579;
t9 = t365 * t41 + t368 * t40 + t607;
t8 = t24 * t362 + t33 * t360;
t7 = t24 * t360 - t33 * t362;
t6 = pkin(2) * t55 - pkin(3) * t114 + pkin(9) * t83 + t26 * t365 + t368 * t66;
t5 = -pkin(9) * t30 + (-pkin(10) * t368 + t490) * t38;
t4 = -pkin(3) * t23 - t512;
t3 = pkin(2) * t16 + pkin(9) * t31 + (-pkin(10) * t365 + t420) * t38;
t2 = -pkin(9) * t23 - t13 * t365 + t15 * t368;
t1 = pkin(2) * t7 - pkin(3) * t33 + pkin(9) * t24 + t13 * t368 + t15 * t365;
t63 = [0, 0, 0, 0, 0, qJDD(1), t387, t388, 0, 0, (t358 * t429 - t440 * t508) * t366, t363 * t340 + (t366 * t316 + (t343 + t383) * t369) * t361, t363 * t314 + (t447 + t369 * (t500 - t415)) * t361, (t358 * t428 + t441 * t508) * t369, t363 * t315 + (t366 * (-t500 + t416) + t445) * t361, t344, (-t292 + pkin(1) * (t337 * t369 + t339 * t366)) * t363 + (t369 * t320 + pkin(1) * t316 + pkin(8) * (t339 * t369 - t447)) * t361, -t320 * t451 - t363 * t293 + pkin(1) * ((t327 * t369 - t338 * t366) * t363 + (t434 * t440 - t429) * t358) + (-t366 * t327 - t445) * t488, pkin(1) * ((-t314 * t369 + t315 * t366) * t363 - (-t359 - t499) * t358 * t435) + (t366 * t314 + t315 * t369) * t488 + t519 * t361, pkin(1) * (t361 * t320 + (-t292 * t369 + t293 * t366) * t363) + t519 * t488, t422 + (t366 * (t303 * t362 - t332 * t453) + t369 * (t303 * t360 + t332 * t452)) * t361, t363 * t301 + (t366 * (-t275 * t362 - t360 * t520) + t369 * (-t275 * t360 + t362 * t520)) * t361, t363 * t280 + (t366 * (-t318 * t360 + t544) + t369 * (t318 * t362 + t545)) * t361, (t330 * t452 + t360 * t380) * t451 + (t330 * t453 - t362 * t380) * t450 - t422, t363 * t375 + (t366 * (t317 * t362 - t460) + t369 * (t317 * t360 + t459)) * t361, t344 + (t366 * (-t330 * t362 + t332 * t360) + t369 * (-t330 * t360 - t332 * t362)) * t361 * t355, (t150 + pkin(1) * (t244 * t369 + t245 * t366)) * t363 + (t366 * (-qJ(3) * t244 - t464) + t369 * (-pkin(2) * t275 + qJ(3) * t245 + t463) - pkin(1) * t275 + pkin(8) * (-t366 * t244 + t245 * t369)) * t361, (t151 + pkin(1) * (t263 * t369 + t264 * t366)) * t363 + (t366 * (-qJ(3) * t263 - t463) + t369 * (-pkin(2) * t520 + qJ(3) * t264 - t464) - pkin(1) * t520 + pkin(8) * (-t366 * t263 + t264 * t369)) * t361, (t492 + pkin(1) * (t229 * t369 + t230 * t366)) * t363 + (t366 * (-qJ(3) * t229 - t126) + t369 * (-pkin(2) * t282 + qJ(3) * t230 + t127) - pkin(1) * t282 + pkin(8) * (-t366 * t229 + t230 * t369)) * t361, (t493 + pkin(1) * (t126 * t369 + t127 * t366)) * t363 + (-qJ(3) * t449 + t369 * (pkin(2) * t270 + qJ(3) * t127) + pkin(1) * t270 + pkin(8) * (t127 * t369 - t449)) * t361, t363 * t221 + (t366 * (t222 * t362 + t424) + t369 * (t222 * t360 - t423)) * t361, t363 * t185 + (t366 * (t187 * t362 + t273 * t360) + t369 * (t187 * t360 - t273 * t362)) * t361, t363 * t213 + (t366 * (t215 * t362 + t227 * t360) + t369 * (t215 * t360 - t227 * t362)) * t361, t363 * t219 + (t366 * (t220 * t362 - t424) + t369 * (t220 * t360 + t423)) * t361, t363 * t214 + (t366 * (t216 * t362 - t223 * t360) + t369 * (t216 * t360 + t223 * t362)) * t361, (t362 * t240 - t360 * t378) * t451 + (t360 * t240 + t362 * t378) * t450 + t363 * t239, (t75 + pkin(1) * (t142 * t369 + t143 * t366)) * t363 + (t366 * (-qJ(3) * t142 + t124 * t362 - t360 * t85) + t369 * (-pkin(2) * t199 + qJ(3) * t143 + t124 * t360 + t362 * t85) - pkin(1) * t199 + pkin(8) * (-t366 * t142 + t143 * t369)) * t361, (t80 + pkin(1) * (t148 * t369 + t149 * t366)) * t363 + (t366 * (-qJ(3) * t148 + t125 * t362 - t360 * t90) + t369 * (-pkin(2) * t208 + qJ(3) * t149 + t125 * t360 + t362 * t90) - pkin(1) * t208 + pkin(8) * (-t366 * t148 + t149 * t369)) * t361, (t46 + pkin(1) * (t136 * t369 + t137 * t366)) * t363 + (t366 * (-qJ(3) * t136 + t362 * t61) + t369 * (qJ(3) * t137 + t360 * t61) + pkin(8) * (-t366 * t136 + t137 * t369) + (t366 * t491 + t369 * t421 - pkin(1)) * t186) * t361, (t37 + pkin(1) * (t366 * t60 + t369 * t59)) * t363 + ((t366 * (-pkin(9) * t362 + t491) + t369 * (-pkin(9) * t360 + t421) - pkin(1)) * t72 + (pkin(8) + qJ(3)) * (-t366 * t59 + t369 * t60)) * t361, t504, -t602, t575, t529, -t605, t528, (t21 + t577) * t363 + (t366 * (t28 * t362 - t360 * t44 - t580) + t369 * (t28 * t360 + t362 * t44 + t578) + t576) * t361, (t22 + t610) * t363 + (t366 * (t29 * t362 - t360 * t45 - t608) + t369 * (t29 * t360 + t362 * t45 - t606) - t609) * t361, (t12 + pkin(1) * (t366 * t58 + t369 * t56)) * t363 + (t366 * (-qJ(3) * t56 + t25 * t362 - t27 * t360) + t369 * (-pkin(2) * t82 + qJ(3) * t58 + t25 * t360 + t27 * t362) - pkin(1) * t82 + pkin(8) * (-t366 * t56 + t369 * t58)) * t361, (t3 + pkin(1) * (t16 * t369 + t17 * t366)) * t363 + (t366 * (-qJ(3) * t16 - t11 * t360 + t362 * t5) + t369 * (-pkin(2) * t30 + qJ(3) * t17 + t11 * t362 + t360 * t5) - pkin(1) * t30 + pkin(8) * (-t366 * t16 + t17 * t369)) * t361, t504, t575, t602, t528, t605, t529, (t10 + t577) * t363 + (t366 * (t19 * t362 - t35 * t360 - t580) + t369 * (t19 * t360 + t35 * t362 + t578) + t576) * t361, (t6 + pkin(1) * (t366 * t57 + t369 * t55)) * t363 + (t366 * (-qJ(3) * t55 + t14 * t362 - t20 * t360) + t369 * (-pkin(2) * t81 + qJ(3) * t57 + t14 * t360 + t20 * t362) - pkin(1) * t81 + pkin(8) * (-t366 * t55 + t369 * t57)) * t361, (t9 - t610) * t363 + (t366 * (t18 * t362 - t32 * t360 + t608) + t369 * (t18 * t360 + t32 * t362 + t606) + t609) * t361, (t1 + pkin(1) * (t366 * t8 + t369 * t7)) * t363 + (t366 * (-qJ(3) * t7 + t2 * t362 - t360 * t4) + t369 * (-pkin(2) * t23 + qJ(3) * t8 + t2 * t360 + t362 * t4) - pkin(1) * t23 + pkin(8) * (-t366 * t7 + t369 * t8)) * t361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t350, t340, t314, t350, t315, t354, -t292, -t293, 0, 0, t302, t301, t280, -t302, t375, t354, t150, t151, t492, t493, t221, t185, t213, t219, t214, t239, t75, t80, t46, t37, t402, -t86, t563, t507, -t104, t509, t21, t22, t12, t3, t402, t563, t86, t509, t104, t507, t10, t6, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, t520, t282, -t270, 0, 0, 0, 0, 0, 0, t199, t208, t186, t72, 0, 0, 0, 0, 0, 0, t549, t98, t82, t30, 0, 0, 0, 0, 0, 0, t549, t81, -t98, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t274, t273, t227, -t274, -t223, -t378, -t122, -t123, 0, 0, t170, t115, t551, t167, -t154, -t527, t386, t385, t384, t404, t170, t551, -t115, -t527, t154, t167, t511, t381, t510, t512; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t461, t523, t516, -t461, -t173, t266, -t64, -t65, 0, 0, t461, t516, -t523, t266, t173, -t461, t374, t396, t304 + t513, t403; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t518, t516, t524, t54;];
tauJ_reg  = t63;
