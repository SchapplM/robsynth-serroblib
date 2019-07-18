% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 05:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRRRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:49:12
% EndTime: 2019-05-08 05:50:27
% DurationCPUTime: 28.60s
% Computational Cost: add. (127570->692), mult. (271914->962), div. (0->0), fcn. (222443->12), ass. (0->453)
t408 = sin(qJ(3));
t412 = cos(qJ(3));
t405 = cos(pkin(6));
t480 = qJD(1) * t405 + qJD(2);
t409 = sin(qJ(2));
t404 = sin(pkin(6));
t508 = qJD(1) * t404;
t492 = t409 * t508;
t378 = -t408 * t492 + t412 * t480;
t379 = t408 * t480 + t412 * t492;
t407 = sin(qJ(4));
t411 = cos(qJ(4));
t357 = -t411 * t378 + t379 * t407;
t354 = qJD(5) + t357;
t566 = t354 ^ 2;
t359 = t378 * t407 + t379 * t411;
t406 = sin(qJ(5));
t410 = cos(qJ(5));
t413 = cos(qJ(2));
t491 = t413 * t508;
t394 = -qJD(3) + t491;
t455 = -qJD(4) + t394;
t334 = t359 * t406 + t410 * t455;
t567 = t334 ^ 2;
t298 = t567 - t566;
t502 = qJD(1) * qJD(2);
t488 = t413 * t502;
t501 = qJDD(1) * t404;
t383 = t404 * t488 + t409 * t501;
t477 = qJDD(1) * t405 + qJDD(2);
t431 = -t412 * t383 - t408 * t477;
t347 = t378 * qJD(3) - t431;
t483 = t383 * t408 - t412 * t477;
t449 = -qJD(3) * t379 - t483;
t486 = t407 * t347 - t411 * t449;
t283 = -t359 * qJD(4) - t486;
t282 = qJDD(5) - t283;
t336 = t410 * t359 - t406 * t455;
t526 = t336 * t334;
t586 = t282 + t526;
t539 = t586 * t406;
t192 = -t298 * t410 + t539;
t305 = t354 * t336;
t284 = -t357 * qJD(4) + t411 * t347 + t407 * t449;
t443 = -qJDD(1) * t413 + t409 * t502;
t439 = t443 * t404;
t428 = qJDD(3) + t439;
t426 = qJDD(4) + t428;
t487 = -t284 * t406 + t410 * t426;
t448 = qJD(5) * t336 - t487;
t218 = -t305 + t448;
t130 = t192 * t407 - t218 * t411;
t134 = t192 * t411 + t218 * t407;
t538 = t586 * t410;
t189 = t298 * t406 + t538;
t85 = t130 * t412 + t134 * t408;
t662 = t404 * (-t413 * t189 + t409 * (t130 * t408 - t134 * t412)) - t405 * t85;
t418 = -t410 * t284 - t406 * t426;
t417 = -t334 * qJD(5) - t418;
t527 = t334 * t354;
t577 = -t527 + t417;
t541 = t577 * t406;
t585 = t305 + t448;
t144 = -t585 * t410 - t541;
t333 = t336 ^ 2;
t582 = t333 - t567;
t109 = t144 * t407 - t411 * t582;
t111 = t144 * t411 + t407 * t582;
t540 = t577 * t410;
t140 = -t585 * t406 + t540;
t70 = t109 * t412 + t111 * t408;
t660 = (t413 * t140 + t409 * (t109 * t408 - t111 * t412)) * t404 - t405 * t70;
t299 = -t333 + t566;
t587 = t282 - t526;
t537 = t587 * t406;
t619 = t299 * t410 + t537;
t576 = t527 + t417;
t536 = t587 * t410;
t618 = -t299 * t406 + t536;
t632 = t407 * t576 + t411 * t618;
t633 = t407 * t618 - t411 * t576;
t640 = t408 * t632 + t412 * t633;
t655 = t405 * t640 + (t409 * (-t408 * t633 + t412 * t632) - t413 * t619) * t404;
t583 = -t333 - t566;
t166 = -t410 * t583 + t539;
t654 = pkin(2) * t166;
t575 = -t566 - t567;
t595 = t410 * t575 - t537;
t614 = t407 * t585 + t411 * t595;
t615 = t407 * t595 - t411 * t585;
t630 = t408 * t614 + t412 * t615;
t653 = pkin(2) * t630;
t652 = pkin(3) * t166;
t651 = pkin(4) * t166;
t650 = pkin(9) * t630;
t649 = pkin(11) * t166;
t168 = t406 * t583 + t538;
t648 = pkin(11) * t168;
t647 = t166 * t413;
t646 = t168 * t407;
t645 = t168 * t411;
t644 = t409 * t166;
t596 = t406 * t575 + t536;
t631 = -t408 * t615 + t412 * t614;
t642 = -pkin(2) * t596 + pkin(9) * t631;
t641 = pkin(1) * (t409 * t631 - t413 * t596);
t639 = pkin(8) * (t409 * t596 + t413 * t631) - pkin(1) * t630;
t637 = pkin(3) * t615;
t636 = pkin(10) * t615;
t634 = -pkin(3) * t596 + pkin(10) * t614;
t627 = pkin(4) * t596;
t626 = pkin(11) * t595;
t625 = pkin(11) * t596;
t624 = -qJ(6) * t406 - pkin(4);
t523 = t354 * t406;
t206 = t334 * t523 - t410 * t448;
t522 = t354 * t410;
t498 = t334 * t522;
t425 = t406 * t448 + t498;
t500 = t407 * t526;
t571 = t411 * t425 - t500;
t499 = t411 * t526;
t572 = t407 * t425 + t499;
t593 = t408 * t571 + t412 * t572;
t617 = t405 * t593 + (t409 * (-t408 * t572 + t412 * t571) - t413 * t206) * t404;
t296 = t336 * t523;
t468 = t296 - t498;
t570 = t282 * t407 + t411 * t468;
t573 = -t411 * t282 + t407 * t468;
t588 = (t334 * t406 + t336 * t410) * t354;
t594 = t408 * t570 + t412 * t573;
t616 = t405 * t594 + (t409 * (-t408 * t573 + t412 * t570) + t413 * t588) * t404;
t581 = t333 + t567;
t613 = pkin(4) * t581;
t612 = qJ(6) * t577;
t609 = t407 * t581;
t316 = t359 * t357;
t584 = -t316 + t426;
t607 = t407 * t584;
t603 = t411 * t581;
t601 = t411 * t584;
t343 = t357 * t455;
t599 = t343 + t284;
t209 = t336 * t522 + t406 * t417;
t210 = t410 * t417 - t296;
t469 = t411 * t210 + t500;
t470 = t407 * t210 - t499;
t517 = t404 * t413;
t518 = t404 * t409;
t569 = t408 * t469 + t412 * t470;
t597 = t405 * t569 + (-t408 * t470 + t412 * t469) * t518 - t209 * t517;
t446 = t455 ^ 2;
t475 = t480 ^ 2;
t561 = sin(qJ(1));
t562 = cos(qJ(1));
t441 = g(1) * t561 - g(2) * t562;
t430 = qJDD(1) * pkin(1) + t441;
t565 = qJD(1) ^ 2;
t505 = t404 * t565;
t421 = pkin(8) * t505 + t430;
t419 = t405 * t421;
t442 = g(1) * t562 + g(2) * t561;
t380 = -pkin(1) * t565 + pkin(8) * t501 - t442;
t484 = -g(3) * t518 + t413 * t380;
t351 = t409 * t419 + t484;
t559 = pkin(2) * t413;
t474 = -pkin(9) * t409 - t559;
t451 = t565 * t474;
t401 = t404 ^ 2;
t510 = t413 * t401;
t579 = -t475 * pkin(2) + t477 * pkin(9);
t415 = t451 * t510 + t351 + t579;
t423 = (-pkin(1) - t559) * qJDD(1) - t441;
t556 = pkin(9) * t413;
t560 = pkin(2) * t409;
t450 = (-t556 + 0.2e1 * t560) * qJD(2);
t554 = t405 * g(3);
t471 = -t383 * pkin(9) - t554;
t473 = -t556 + t560;
t555 = t404 * pkin(8);
t416 = (((t405 * t473 - t555) * qJD(1) + t450) * qJD(1) + t423) * t404 + t471;
t267 = t408 * t416 + t412 * t415;
t362 = -pkin(3) * t394 - pkin(10) * t379;
t376 = t378 ^ 2;
t229 = -t376 * pkin(3) + pkin(10) * t449 + t394 * t362 + t267;
t368 = t378 * t394;
t326 = t368 + t347;
t521 = t378 * t379;
t422 = t428 + t521;
t414 = pkin(3) * t422 - pkin(10) * t326 - t408 * t415 + t412 * t416;
t152 = t411 * t229 + t407 * t414;
t311 = pkin(4) * t357 - pkin(11) * t359;
t125 = -pkin(4) * t446 + pkin(11) * t426 - t357 * t311 + t152;
t485 = t409 * t380 - t413 * t419;
t319 = -t477 * pkin(2) - t475 * pkin(9) + (t413 * g(3) + t451 * t518) * t404 + t485;
t251 = -t449 * pkin(3) - t376 * pkin(10) + t379 * t362 + t319;
t440 = t455 * t359;
t149 = -t599 * pkin(11) + (-t283 - t440) * pkin(4) + t251;
t90 = t125 * t406 - t410 * t149;
t91 = t410 * t125 + t406 * t149;
t44 = t406 * t90 + t410 * t91;
t590 = t408 * t422;
t589 = t412 * t422;
t507 = qJD(6) * t354;
t349 = 0.2e1 * t507;
t286 = pkin(5) * t334 - qJ(6) * t336;
t479 = -t282 * qJ(6) + t334 * t286 - t91;
t465 = t349 - t479;
t61 = -pkin(5) * t566 + t465;
t62 = -t282 * pkin(5) - qJ(6) * t566 + t286 * t336 + qJDD(6) + t90;
t38 = t406 * t62 + t410 * t61;
t151 = t229 * t407 - t411 * t414;
t124 = -t426 * pkin(4) - t446 * pkin(11) + t311 * t359 + t151;
t427 = t448 * pkin(5) + t124 - t612;
t82 = (pkin(5) * t354 - 0.2e1 * qJD(6)) * t336 + t427;
t27 = t38 * t407 - t411 * t82;
t433 = pkin(11) * t38 + (-pkin(5) * t410 + t624) * t82;
t580 = pkin(3) * t27 + t433;
t350 = g(3) * t517 + t485;
t578 = t409 * t350 + t413 * t351;
t574 = -pkin(5) * (t583 + t566) + qJ(6) * t586 - t479;
t322 = (qJD(3) + t394) * t379 + t483;
t355 = t357 ^ 2;
t356 = t359 ^ 2;
t377 = t379 ^ 2;
t391 = t394 ^ 2;
t94 = -t151 * t411 + t152 * t407;
t563 = pkin(3) * t94;
t253 = t359 * t394 + t486;
t256 = -t343 + t284;
t182 = -t253 * t407 - t256 * t411;
t558 = pkin(3) * t182;
t557 = pkin(4) * t407;
t552 = t408 * t94;
t551 = t412 * t94;
t550 = -pkin(4) * t124 + pkin(11) * t44;
t548 = qJ(6) * t410;
t543 = t576 * t406;
t542 = t576 * t410;
t535 = t251 * t407;
t534 = t251 * t411;
t301 = -t316 - t426;
t530 = t301 * t407;
t529 = t301 * t411;
t528 = t319 * t412;
t340 = -t428 + t521;
t525 = t340 * t408;
t524 = t340 * t412;
t520 = t394 * t408;
t519 = t394 * t412;
t516 = t405 * t409;
t120 = t406 * t124;
t515 = t408 * t319;
t506 = t401 * t565;
t393 = t413 * t409 * t506;
t382 = t393 + t477;
t513 = t409 * t382;
t121 = t410 * t124;
t381 = -t393 + t477;
t511 = t413 * t381;
t504 = -qJD(3) + t394;
t497 = t413 * t316;
t496 = t413 * t521;
t225 = (qJD(5) + t354) * t334 + t418;
t495 = pkin(4) * t225 + t120 - t648;
t494 = -pkin(4) * t585 - t121 + t626;
t493 = -pkin(4) * t411 - pkin(3);
t402 = t409 ^ 2;
t490 = t402 * t506;
t403 = t413 ^ 2;
t489 = t403 * t506;
t95 = t151 * t407 + t411 * t152;
t266 = t408 * (t430 * t516 + t484 + t579) - t412 * t471 + (-t412 * t423 + (((t408 * t413 * t474 + t412 * pkin(8)) * t404 + (t408 * t409 * pkin(8) - t412 * t473) * t405) * qJD(1) - t412 * t450) * qJD(1)) * t404;
t201 = t266 * t408 + t412 * t267;
t143 = -t218 * t410 + t543;
t52 = (t581 - t566) * pkin(5) + t465;
t54 = qJ(6) * t581 + t62;
t482 = pkin(11) * t143 + t406 * t54 + t410 * t52 + t613;
t220 = (-qJD(5) + t354) * t336 + t487;
t145 = t220 * t410 + t543;
t481 = pkin(11) * t145 + t44 + t613;
t39 = -t124 * t411 + t407 * t44;
t478 = pkin(3) * t39 + t550;
t328 = -t356 - t446;
t264 = t328 * t411 + t530;
t476 = pkin(3) * t264 - t152;
t472 = -pkin(5) * t62 + qJ(6) * t61;
t467 = -pkin(5) * t576 - qJ(6) * t218;
t43 = t406 * t91 - t410 * t90;
t463 = t266 * t412 - t267 * t408;
t461 = t401 * t409 * t488;
t460 = -pkin(1) + t474;
t459 = t494 + t637;
t116 = t225 * t411 - t646;
t458 = pkin(3) * t116 + t495;
t308 = -t446 - t355;
t245 = t308 * t407 + t601;
t456 = pkin(3) * t245 - t151;
t101 = t143 * t407 + t603;
t454 = pkin(3) * t101 + t482;
t102 = t145 * t407 + t603;
t453 = pkin(3) * t102 + t481;
t424 = 0.2e1 * qJD(6) * t336 - t427;
t57 = -pkin(5) * t305 + t424 + t612;
t452 = pkin(4) * t577 + pkin(5) * t540 + t406 * t57 + t648;
t58 = (-t585 - t305) * pkin(5) + t424;
t447 = t410 * t58 + t624 * t585 + t626;
t445 = t480 * t508;
t112 = t411 * t577 + t646;
t438 = pkin(3) * t112 + t452;
t437 = t407 * t343;
t436 = t407 * t440;
t435 = t411 * t343;
t434 = t411 * t440;
t432 = t447 + t637;
t420 = pkin(5) * t587 + qJ(6) * t575 - t62;
t387 = t413 * t445;
t386 = t409 * t445;
t385 = (t402 - t403) * t506;
t384 = -t475 - t489;
t373 = -t490 - t475;
t369 = t404 * t421 + t554;
t367 = -t386 - t439;
t366 = t386 - t439;
t365 = -t387 + t383;
t364 = -t377 + t391;
t363 = t376 - t391;
t361 = -t377 - t391;
t360 = t377 - t376;
t348 = -t391 - t376;
t339 = -t356 + t446;
t338 = t355 - t446;
t337 = t376 + t377;
t330 = (-t378 * t408 + t379 * t412) * t394;
t327 = t378 * t504 + t431;
t325 = -t368 + t347;
t323 = t379 * t504 - t483;
t321 = t347 * t408 - t379 * t519;
t320 = t378 * t520 + t412 * t449;
t315 = t356 - t355;
t314 = t363 * t408 - t524;
t313 = t364 * t412 + t590;
t310 = -t361 * t408 + t524;
t309 = t361 * t412 + t525;
t295 = t348 * t412 - t590;
t294 = t348 * t408 + t589;
t293 = t435 - t436;
t292 = t437 + t434;
t285 = -t355 - t356;
t276 = -t322 * t412 + t326 * t408;
t274 = t323 * t408 + t325 * t412;
t271 = t338 * t411 + t530;
t270 = -t339 * t407 + t601;
t269 = t338 * t407 - t529;
t268 = t339 * t411 + t607;
t265 = -t328 * t407 + t529;
t252 = (0.2e1 * qJD(4) - t394) * t359 + t486;
t250 = t411 * t284 + t436;
t249 = t407 * t284 - t434;
t248 = -t407 * t283 - t435;
t247 = t411 * t283 - t437;
t246 = t308 * t411 - t607;
t240 = t292 * t412 + t293 * t408;
t236 = pkin(2) * t327 + pkin(9) * t310 + t515;
t230 = pkin(2) * t323 + pkin(9) * t295 - t528;
t216 = t269 * t412 + t271 * t408;
t215 = t268 * t412 + t270 * t408;
t199 = -t264 * t408 + t265 * t412;
t198 = t264 * t412 + t265 * t408;
t185 = -pkin(10) * t264 + t534;
t184 = -t253 * t411 + t256 * t407;
t183 = -t252 * t411 - t407 * t599;
t181 = -t252 * t407 + t411 * t599;
t180 = -pkin(10) * t245 + t535;
t179 = t249 * t412 + t250 * t408;
t178 = t247 * t412 + t248 * t408;
t173 = -t245 * t408 + t246 * t412;
t172 = t245 * t412 + t246 * t408;
t165 = -pkin(2) * t319 + pkin(9) * t201;
t154 = pkin(2) * t337 + pkin(9) * t276 + t201;
t153 = -pkin(3) * t599 + pkin(10) * t265 + t535;
t146 = -pkin(3) * t252 + pkin(10) * t246 - t534;
t141 = t220 * t406 - t542;
t139 = -t218 * t406 - t542;
t118 = -t225 * t407 - t645;
t114 = -t407 * t577 + t645;
t107 = -t182 * t408 + t184 * t412;
t106 = t182 * t412 + t184 * t408;
t105 = t181 * t412 + t183 * t408;
t104 = t145 * t411 - t609;
t103 = t143 * t411 - t609;
t96 = t121 + t649;
t93 = t120 - t625;
t92 = -pkin(4) * t139 - t467;
t87 = -pkin(3) * t251 + pkin(10) * t95;
t80 = -pkin(2) * t599 + pkin(9) * t199 + t153 * t412 + t185 * t408;
t78 = -t116 * t408 + t118 * t412;
t76 = t116 * t412 + t118 * t408;
t74 = -t112 * t408 + t114 * t412;
t72 = t112 * t412 + t114 * t408;
t71 = -pkin(10) * t182 - t94;
t68 = -pkin(2) * t252 + pkin(9) * t173 + t146 * t412 + t180 * t408;
t67 = -t102 * t408 + t104 * t412;
t66 = -t101 * t408 + t103 * t412;
t65 = t102 * t412 + t104 * t408;
t64 = t101 * t412 + t103 * t408;
t63 = -pkin(3) * t285 + pkin(10) * t184 + t95;
t60 = t91 + t651;
t59 = t90 - t627;
t50 = t412 * t95 - t552;
t49 = t408 * t95 + t551;
t48 = -t420 - t627;
t47 = -t406 * t58 - t548 * t585 - t625;
t46 = -pkin(5) * t541 + t410 * t57 - t649;
t45 = -0.2e1 * t507 - t574 - t651;
t41 = -pkin(11) * t141 - t43;
t40 = t124 * t407 + t411 * t44;
t37 = t406 * t61 - t410 * t62;
t35 = -pkin(10) * t116 - t407 * t60 + t411 * t96;
t34 = -t407 * t59 + t411 * t93 - t636;
t33 = pkin(10) * t118 + t407 * t96 + t411 * t60 + t652;
t32 = t407 * t93 + t411 * t59 + t634;
t31 = -pkin(2) * t285 + pkin(9) * t107 + t408 * t71 + t412 * t63;
t30 = -pkin(11) * t139 - t406 * t52 + t410 * t54;
t29 = -pkin(10) * t102 + t141 * t557 + t41 * t411;
t28 = t38 * t411 + t407 * t82;
t26 = -pkin(2) * t251 + pkin(9) * t50 - pkin(10) * t552 + t412 * t87;
t25 = pkin(10) * t104 + t141 * t493 + t407 * t41;
t24 = -t407 * t48 + t411 * t47 - t636;
t23 = -pkin(10) * t112 - t407 * t45 + t411 * t46;
t22 = t407 * t47 + t411 * t48 + t634;
t21 = pkin(10) * t114 + t407 * t46 + t411 * t45 - t652;
t20 = -pkin(11) * t37 + (pkin(5) * t406 - t548) * t82;
t19 = -pkin(10) * t101 + t30 * t411 - t407 * t92;
t18 = -t39 * t408 + t40 * t412;
t17 = t39 * t412 + t40 * t408;
t16 = -pkin(4) * t37 - t472;
t15 = -pkin(3) * t139 + pkin(10) * t103 + t30 * t407 + t411 * t92;
t14 = -pkin(10) * t39 + (-pkin(11) * t411 + t557) * t43;
t13 = pkin(9) * t78 + t33 * t412 + t35 * t408 + t654;
t12 = t32 * t412 + t34 * t408 + t642;
t11 = -t27 * t408 + t28 * t412;
t10 = t27 * t412 + t28 * t408;
t9 = pkin(10) * t40 + (-pkin(11) * t407 + t493) * t43;
t8 = -pkin(2) * t141 + pkin(9) * t67 + t25 * t412 + t29 * t408;
t7 = t22 * t412 + t24 * t408 + t642;
t6 = pkin(9) * t74 + t21 * t412 + t23 * t408 - t654;
t5 = -pkin(2) * t139 + pkin(9) * t66 + t15 * t412 + t19 * t408;
t4 = -pkin(10) * t27 - t16 * t407 + t20 * t411;
t3 = -pkin(3) * t37 + pkin(10) * t28 + t16 * t411 + t20 * t407;
t2 = -pkin(2) * t43 + pkin(9) * t18 + t14 * t408 + t412 * t9;
t1 = -pkin(2) * t37 + pkin(9) * t11 + t3 * t412 + t4 * t408;
t36 = [0, 0, 0, 0, 0, qJDD(1), t441, t442, 0, 0, t383 * t518 + t461, t405 * t385 + (t409 * t367 + t413 * (t387 + t383)) * t404, t405 * t365 + (t513 + t413 * (t475 - t490)) * t404, -t443 * t510 - t461, t405 * t366 + (t409 * (-t475 + t489) + t511) * t404, t405 * t477, (-t350 + pkin(1) * (t382 * t413 + t384 * t409)) * t405 + (t413 * t369 + pkin(1) * t367 + pkin(8) * (t384 * t413 - t513)) * t404, -t369 * t518 - t405 * t351 + pkin(1) * (-t381 * t516 + t405 * t413 * t373 - t404 * (t480 * t491 + t383)) + (-t409 * t373 - t511) * t555, pkin(1) * ((-t365 * t413 + t366 * t409) * t405 - (-t402 - t403) * t401 * t505) + (t409 * t365 + t366 * t413) * t555 + t578 * t404, pkin(1) * (t404 * t369 + (-t350 * t413 + t351 * t409) * t405) + t578 * t555, t405 * t321 + (t409 * (t347 * t412 + t379 * t520) + t496) * t404, t405 * t274 + (t409 * (t323 * t412 - t325 * t408) - t413 * t360) * t404, t405 * t313 + (t409 * (-t364 * t408 + t589) - t413 * t326) * t404, t405 * t320 + (t409 * (t378 * t519 - t408 * t449) - t496) * t404, t405 * t314 + (t409 * (t363 * t412 + t525) + t413 * t322) * t404, -t428 * t517 + t405 * t330 + (-t378 * t412 - t379 * t408) * t394 * t518, (-pkin(9) * t294 + t515) * t518 + (-pkin(2) * t294 + t266) * t517 + t405 * t230 + pkin(1) * (-t404 * t294 + (t295 * t409 + t323 * t413) * t405) + (t295 * t413 - t409 * t323) * t555, (t236 + pkin(1) * (t310 * t409 + t327 * t413)) * t405 + (t409 * (-pkin(9) * t309 + t528) + t413 * (-pkin(2) * t309 + t267) - pkin(1) * t309 + pkin(8) * (t310 * t413 - t409 * t327)) * t404, (t154 + pkin(1) * (t276 * t409 + t337 * t413)) * t405 + (t409 * t463 + pkin(8) * (t276 * t413 - t409 * t337) + t460 * (-t322 * t408 - t326 * t412)) * t404, (t165 + pkin(1) * (t201 * t409 - t319 * t413)) * t405 + (pkin(8) * (t201 * t413 + t409 * t319) - t460 * t463) * t404, t405 * t179 + (t409 * (-t249 * t408 + t250 * t412) - t497) * t404, t405 * t105 + (t409 * (-t181 * t408 + t183 * t412) - t413 * t315) * t404, t405 * t215 + (t409 * (-t268 * t408 + t270 * t412) - t413 * t256) * t404, t405 * t178 + (t409 * (-t247 * t408 + t248 * t412) + t497) * t404, t405 * t216 + (t409 * (-t269 * t408 + t271 * t412) + t413 * t253) * t404, t405 * t240 + (t409 * (-t292 * t408 + t293 * t412) - t426 * t413) * t404, (t68 + pkin(1) * (t173 * t409 - t252 * t413)) * t405 + (t409 * (-pkin(9) * t172 - t146 * t408 + t180 * t412) + t413 * (-pkin(2) * t172 - t456) - pkin(1) * t172 + pkin(8) * (t173 * t413 + t409 * t252)) * t404, (t80 + pkin(1) * (t199 * t409 - t413 * t599)) * t405 + (t409 * (-pkin(9) * t198 - t153 * t408 + t185 * t412) + t413 * (-pkin(2) * t198 - t476) - pkin(1) * t198 + pkin(8) * (t199 * t413 + t409 * t599)) * t404, (t31 + pkin(1) * (t107 * t409 - t285 * t413)) * t405 + (t409 * (-pkin(9) * t106 - t408 * t63 + t412 * t71) + t413 * (-pkin(2) * t106 - t558) - pkin(1) * t106 + pkin(8) * (t107 * t413 + t409 * t285)) * t404, (t26 + pkin(1) * (-t251 * t413 + t409 * t50)) * t405 + (t409 * (-pkin(9) * t49 - pkin(10) * t551 - t408 * t87) + t413 * (-pkin(2) * t49 - t563) - pkin(1) * t49 + pkin(8) * (t409 * t251 + t413 * t50)) * t404, t597, -t660, t655, t617, t662, t616, (t12 + t641) * t405 + (t409 * (-t32 * t408 + t34 * t412 - t650) + t413 * (-t459 - t653) + t639) * t404, (t13 + pkin(1) * (t409 * t78 + t647)) * t405 + (t409 * (-pkin(9) * t76 - t33 * t408 + t35 * t412) + t413 * (-pkin(2) * t76 - t458) - pkin(1) * t76 + pkin(8) * (t413 * t78 - t644)) * t404, (t8 + pkin(1) * (-t141 * t413 + t409 * t67)) * t405 + (t409 * (-pkin(9) * t65 - t25 * t408 + t29 * t412) + t413 * (-pkin(2) * t65 - t453) - pkin(1) * t65 + pkin(8) * (t409 * t141 + t413 * t67)) * t404, (t2 + pkin(1) * (t18 * t409 - t413 * t43)) * t405 + (t409 * (-pkin(9) * t17 + t14 * t412 - t408 * t9) + t413 * (-pkin(2) * t17 - t478) - pkin(1) * t17 + pkin(8) * (t18 * t413 + t409 * t43)) * t404, t597, t655, t660, t616, -t662, t617, (t7 + t641) * t405 + (t409 * (-t22 * t408 + t24 * t412 - t650) + t413 * (-t432 - t653) + t639) * t404, (t5 + pkin(1) * (-t139 * t413 + t409 * t66)) * t405 + (t409 * (-pkin(9) * t64 - t15 * t408 + t19 * t412) + t413 * (-pkin(2) * t64 - t454) - pkin(1) * t64 + pkin(8) * (t409 * t139 + t413 * t66)) * t404, (t6 + pkin(1) * (t409 * t74 - t647)) * t405 + (t409 * (-pkin(9) * t72 - t21 * t408 + t23 * t412) + t413 * (-pkin(2) * t72 - t438) - pkin(1) * t72 + pkin(8) * (t413 * t74 + t644)) * t404, (t1 + pkin(1) * (t11 * t409 - t37 * t413)) * t405 + (t409 * (-pkin(9) * t10 - t3 * t408 + t4 * t412) + t413 * (-pkin(2) * t10 - t580) - pkin(1) * t10 + pkin(8) * (t11 * t413 + t409 * t37)) * t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t393, t385, t365, t393, t366, t477, -t350, -t351, 0, 0, t321, t274, t313, t320, t314, t330, t230, t236, t154, t165, t179, t105, t215, t178, t216, t240, t68, t80, t31, t26, t569, t70, t640, t593, -t85, t594, t12, t13, t8, t2, t569, t640, -t70, t594, t85, t593, t7, t5, t6, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t521, t360, t326, t521, -t322, t428, -t266, -t267, 0, 0, t316, t315, t256, -t316, -t253, t426, t456, t476, t558, t563, t209, t140, t619, t206, t189, -t588, t459, t458, t453, t478, t209, t619, -t140, -t588, -t189, t206, t432, t454, t438, t580; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, t315, t256, -t316, -t253, t426, -t151, -t152, 0, 0, t209, t140, t619, t206, t189, -t588, t494, t495, t481, t550, t209, t619, -t140, -t588, -t189, t206, t447, t482, t452, t433; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t526, t582, t576, -t526, -t218, t282, -t90, -t91, 0, 0, t526, t576, -t582, t282, t218, -t526, t420, t467, t349 + t574, t472; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t587, t576, t583, t62;];
tauJ_reg  = t36;