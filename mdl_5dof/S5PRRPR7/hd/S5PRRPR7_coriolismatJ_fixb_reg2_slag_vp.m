% Calculate inertial parameters regressor of coriolis matrix for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPR7_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:41
% EndTime: 2019-12-05 16:38:00
% DurationCPUTime: 7.42s
% Computational Cost: add. (5720->491), mult. (14544->816), div. (0->0), fcn. (15707->10), ass. (0->360)
t357 = sin(qJ(3));
t351 = t357 ^ 2;
t360 = cos(qJ(3));
t352 = t360 ^ 2;
t338 = t352 - t351;
t353 = sin(pkin(10));
t355 = cos(pkin(10));
t359 = cos(qJ(5));
t512 = t359 * t360;
t356 = sin(qJ(5));
t516 = t356 * t357;
t295 = t355 * t516 + t512;
t560 = -t359 / 0.2e1;
t435 = t295 * t560;
t347 = t359 * t357;
t515 = t356 * t360;
t298 = t355 * t347 - t515;
t534 = t298 * t356;
t377 = (t435 - t534 / 0.2e1) * t353;
t581 = t377 * qJD(5);
t349 = t353 ^ 2;
t529 = t349 * t359;
t476 = t356 * t529;
t580 = qJD(2) * t377 - qJD(3) * t476;
t503 = qJD(2) * t298;
t579 = -qJD(3) * t377 + t295 * t503;
t354 = sin(pkin(5));
t361 = cos(qJ(2));
t525 = t354 * t361;
t472 = t353 * t525;
t555 = cos(pkin(5));
t431 = t555 * t357;
t358 = sin(qJ(2));
t526 = t354 * t358;
t473 = t360 * t526;
t294 = t431 + t473;
t537 = t294 * t355;
t227 = -t472 + t537;
t293 = t357 * t526 - t555 * t360;
t540 = t293 * t356;
t141 = t227 * t359 + t540;
t578 = -t141 / 0.2e1;
t520 = t355 * t361;
t471 = t354 * t520;
t538 = t294 * t353;
t226 = t471 + t538;
t577 = t226 / 0.2e1;
t267 = -t355 * t526 + t360 * t472;
t576 = -t267 / 0.2e1;
t575 = -t293 / 0.2e1;
t574 = -t295 / 0.2e1;
t573 = t295 / 0.2e1;
t572 = -t298 / 0.2e1;
t571 = t298 / 0.2e1;
t299 = t355 * t512 + t516;
t570 = t299 / 0.2e1;
t569 = t347 / 0.2e1;
t568 = -t353 / 0.2e1;
t567 = t353 / 0.2e1;
t566 = -t355 / 0.2e1;
t565 = t355 / 0.2e1;
t564 = -t356 / 0.2e1;
t563 = t356 / 0.2e1;
t562 = -t357 / 0.2e1;
t561 = t357 / 0.2e1;
t559 = t359 / 0.2e1;
t558 = -t360 / 0.2e1;
t557 = t294 * pkin(7);
t556 = t360 * pkin(7);
t404 = pkin(4) * t353 - pkin(8) * t355 + pkin(7);
t291 = t404 * t360;
t513 = t359 * t291;
t332 = t357 * pkin(3) - t360 * qJ(4);
t312 = t353 * t332;
t244 = t312 + (-pkin(7) * t355 + pkin(8)) * t357;
t519 = t356 * t244;
t146 = t513 - t519;
t554 = t146 * t359;
t514 = t359 * t244;
t517 = t356 * t291;
t147 = t514 + t517;
t553 = t147 * t356;
t523 = t355 * t356;
t535 = t294 * t359;
t184 = t293 * t523 + t535;
t552 = t184 * t359;
t521 = t355 * t359;
t536 = t294 * t356;
t185 = -t293 * t521 + t536;
t551 = t185 * t356;
t470 = t357 * t525;
t268 = (t353 * t358 + t360 * t520) * t354;
t544 = t268 * t356;
t196 = t359 * t470 - t544;
t550 = t196 * t355;
t543 = t268 * t359;
t197 = t356 * t470 + t543;
t549 = t197 * t355;
t548 = t226 * t267;
t547 = t226 * t353;
t546 = t227 * t355;
t405 = -pkin(4) * t355 - pkin(8) * t353 - pkin(3);
t266 = qJ(4) * t521 + t356 * t405;
t545 = t266 * t359;
t422 = -t360 * pkin(3) - t357 * qJ(4);
t319 = -pkin(2) + t422;
t343 = t355 * t556;
t270 = t353 * t319 + t343;
t542 = t270 * t355;
t541 = t293 * t353;
t539 = t293 * t359;
t140 = t227 * t356 - t539;
t30 = -t140 * t184 + t141 * t185 - t226 * t541;
t533 = t30 * qJD(1);
t31 = -t140 * t196 + t141 * t197 + t548;
t532 = t31 * qJD(1);
t531 = t349 * qJ(4);
t530 = t349 * t356;
t528 = t353 * t357;
t527 = t353 * t360;
t524 = t355 * t332;
t522 = t355 * t357;
t265 = qJ(4) * t523 - t359 * t405;
t518 = t356 * t265;
t49 = (t294 - t546 - t547) * t293;
t511 = t49 * qJD(1);
t50 = t227 * t268 + t293 * t470 + t548;
t510 = t50 * qJD(1);
t376 = (t435 + t534 / 0.2e1) * t355;
t297 = t355 * t515 - t347;
t386 = t297 * t564 + t299 * t560;
t79 = t376 - t386;
t509 = t79 * qJD(2);
t99 = -t295 * t299 - t297 * t298;
t508 = t99 * qJD(2);
t350 = t355 ^ 2;
t335 = t350 + t349;
t474 = t353 * t347;
t475 = t353 * t516;
t150 = -t295 * t474 + t298 * t475;
t507 = qJD(2) * t150;
t219 = t295 * t522 + t351 * t530;
t506 = qJD(2) * t219;
t220 = t298 * t522 + t351 * t529;
t505 = qJD(2) * t220;
t504 = qJD(2) * t295;
t502 = qJD(2) * t357;
t501 = qJD(2) * t358;
t500 = qJD(2) * t360;
t499 = qJD(3) * t353;
t498 = qJD(4) * t355;
t497 = qJD(4) * t360;
t496 = qJD(5) * t355;
t122 = (t293 * t357 + t294 * t360 - t526) * t525;
t495 = t122 * qJD(1);
t433 = -t512 / 0.2e1;
t162 = (-t349 / 0.2e1 - 0.1e1 / 0.2e1) * t516 + (t573 + t433) * t355;
t494 = t162 * qJD(2);
t436 = -t515 / 0.2e1;
t163 = t569 + t349 * t569 + (t436 + t572) * t355;
t493 = t163 * qJD(2);
t169 = (-t295 * t360 - t297 * t357) * t353;
t492 = t169 * qJD(2);
t170 = (t298 * t360 + t299 * t357) * t353;
t491 = t170 * qJD(2);
t490 = t293 * qJD(3);
t489 = t294 * qJD(3);
t302 = t335 * t351;
t488 = t302 * qJD(2);
t305 = t338 * t353;
t487 = t305 * qJD(2);
t306 = t338 * t355;
t486 = t306 * qJD(2);
t485 = t335 * qJD(3);
t484 = t338 * qJD(2);
t483 = t355 * qJD(3);
t482 = t357 * qJD(3);
t481 = t360 * qJD(3);
t480 = pkin(7) * t528;
t479 = pkin(2) * t502;
t478 = pkin(2) * t500;
t477 = t556 / 0.2e1;
t468 = t353 * t502;
t467 = t355 * t502;
t466 = t353 * t483;
t465 = t356 * t499;
t464 = t359 * t499;
t463 = qJD(4) * t528;
t462 = t357 * t497;
t461 = qJD(5) * t528;
t460 = t356 * t496;
t459 = t359 * t496;
t458 = qJD(2) * t525;
t457 = t357 * t481;
t341 = t357 * t500;
t456 = t226 * t574;
t455 = t297 * t577;
t454 = t226 * t571;
t453 = t226 * t570;
t452 = -t547 / 0.2e1;
t451 = t547 / 0.2e1;
t450 = t226 * t565;
t449 = t267 * t567;
t342 = pkin(7) * t527;
t269 = t319 * t355 - t342;
t241 = pkin(4) * t360 - t269;
t448 = t241 * t575;
t447 = t540 / 0.2e1;
t446 = t539 / 0.2e1;
t445 = -t528 / 0.2e1;
t444 = t528 / 0.2e1;
t443 = t353 * t559;
t442 = -t527 / 0.2e1;
t441 = t527 / 0.2e1;
t440 = -t525 / 0.2e1;
t439 = t525 / 0.2e1;
t438 = t522 / 0.2e1;
t437 = -t516 / 0.2e1;
t434 = -t347 / 0.2e1;
t432 = t512 / 0.2e1;
t428 = t360 * t466;
t427 = qJD(5) * t476;
t426 = t355 * t341;
t425 = t357 * t440;
t424 = t357 * t439;
t423 = t294 / 0.2e1 + t452;
t421 = -qJD(5) - t468;
t242 = -pkin(8) * t360 + t270;
t382 = t404 * t357;
t142 = t356 * t242 - t359 * t382;
t143 = t359 * t242 + t356 * t382;
t243 = -t524 + (-pkin(7) * t353 - pkin(4)) * t357;
t363 = -t140 * t146 / 0.2e1 + t141 * t147 / 0.2e1 - t184 * t142 / 0.2e1 + t185 * t143 / 0.2e1 + t243 * t577;
t395 = t196 * t265 / 0.2e1 - t197 * t266 / 0.2e1;
t2 = (qJ(4) * t576 + t448) * t353 + t363 + t395;
t21 = -t142 * t146 + t143 * t147 + t241 * t243;
t420 = t2 * qJD(1) + t21 * qJD(2);
t20 = t142 * t299 - t143 * t297 - t146 * t298 - t147 * t295;
t367 = t140 * t570 + t184 * t572 + t185 * t574 + t297 * t578;
t394 = t196 * t560 + t197 * t564;
t379 = t394 * t353;
t5 = t379 - t367;
t419 = -t5 * qJD(1) + t20 * qJD(2);
t389 = t269 * t568 + t542 / 0.2e1;
t123 = t477 - t389;
t390 = t268 * t565 + t449;
t378 = t390 * qJ(4);
t276 = t480 + t524;
t277 = -pkin(7) * t522 + t312;
t393 = t276 * t577 - t227 * t277 / 0.2e1;
t13 = pkin(3) * t425 - t123 * t293 + t557 * t562 + t378 + t393;
t74 = pkin(7) ^ 2 * t357 * t360 + t269 * t276 + t270 * t277;
t418 = -t13 * qJD(1) + t74 * qJD(2);
t371 = t140 * t558 + t184 * t561 + t293 * t574;
t22 = t455 + t550 / 0.2e1 + (t267 * t564 + t371) * t353;
t34 = t241 * t297 + t243 * t295 + (-t142 * t360 + t146 * t357) * t353;
t417 = t22 * qJD(1) + t34 * qJD(2);
t370 = t141 * t558 + t185 * t562 + t293 * t572;
t25 = t453 - t549 / 0.2e1 + (t267 * t560 + t370) * t353;
t35 = -t241 * t299 - t243 * t298 + (t143 * t360 + t147 * t357) * t353;
t416 = t25 * qJD(1) - t35 * qJD(2);
t398 = t140 * t563 + t141 * t559;
t366 = (-t398 * t353 + t450) * t357;
t29 = t366 + t394;
t40 = t142 * t475 + t143 * t474 - t241 * t522;
t415 = qJD(1) * t29 - qJD(2) * t40;
t41 = t456 + t543 / 0.2e1 + (t140 * t567 + t356 * t439) * t357;
t70 = -t142 * t528 + t241 * t295;
t414 = qJD(1) * t41 - qJD(2) * t70;
t43 = t454 + t544 / 0.2e1 + (t141 * t568 + t359 * t440) * t357;
t71 = -t143 * t528 + t241 * t298;
t413 = qJD(1) * t43 + qJD(2) * t71;
t52 = (t360 * t577 - t268 / 0.2e1) * t355 + (t227 * t558 + t576) * t353;
t63 = (t269 * t360 + t276 * t357) * t355 + (t270 * t360 + t277 * t357) * t353;
t412 = -t52 * qJD(1) + t63 * qJD(2);
t411 = -t276 * t353 + t277 * t355;
t118 = -t269 * t357 + (t276 - 0.2e1 * t480) * t360;
t58 = (t355 * t440 + t577 - t538 / 0.2e1) * t357;
t410 = t58 * qJD(1) + t118 * qJD(2);
t119 = t277 * t360 + (-t270 + 0.2e1 * t343) * t357;
t57 = (t353 * t439 + t227 / 0.2e1 - t537 / 0.2e1) * t357;
t409 = t57 * qJD(1) - t119 * qJD(2);
t134 = (t269 * t355 + t270 * t353) * t357;
t391 = t227 * t568 + t450;
t80 = (t439 - t391) * t357;
t408 = qJD(1) * t80 + qJD(2) * t134;
t151 = t295 ^ 2 - t298 ^ 2;
t97 = (-t295 * t356 + t298 * t359) * t353;
t407 = qJD(2) * t151 - qJD(3) * t97;
t303 = (t356 ^ 2 - t359 ^ 2) * t349;
t406 = qJD(2) * t97 - qJD(3) * t303;
t403 = t353 * (-qJD(5) + t483);
t400 = t355 * t437 + t573;
t173 = (t433 + t400) * t353;
t304 = t335 * t356;
t402 = qJD(2) * t173 + qJD(3) * t304;
t399 = t355 * t434 + t571;
t174 = (t515 / 0.2e1 + t399) * t353;
t307 = t335 * t359;
t401 = qJD(2) * t174 + qJD(3) * t307;
t397 = -t554 / 0.2e1 - t553 / 0.2e1;
t396 = -t552 / 0.2e1 - t551 / 0.2e1;
t392 = t451 + t546 / 0.2e1;
t388 = -t519 / 0.2e1 + t513 / 0.2e1;
t387 = -t517 / 0.2e1 - t514 / 0.2e1;
t383 = t468 - t483;
t381 = -qJD(5) - t383;
t102 = t531 + (t518 + t545) * t355;
t368 = t398 * t355 + t451;
t27 = t368 + t396;
t362 = (t241 / 0.2e1 + (-t545 / 0.2e1 - t518 / 0.2e1 + qJ(4) * t565) * t357) * t353 + (t142 * t563 + t143 * t559) * t355;
t8 = t362 + t397;
t380 = -qJD(1) * t27 - qJD(2) * t8 - qJD(3) * t102;
t212 = -qJ(4) * t530 - t265 * t355;
t365 = (qJ(4) * t574 + t241 * t564 + t265 * t561) * t353 + t142 * t566;
t36 = t365 - t387;
t46 = -t423 * t356 + (t446 + t140 / 0.2e1) * t355;
t375 = -qJD(1) * t46 + qJD(2) * t36 + qJD(3) * t212;
t213 = -qJ(4) * t529 - t266 * t355;
t364 = (qJ(4) * t571 + t241 * t559 + t266 * t562) * t353 + t143 * t565;
t38 = t364 - t388;
t45 = t423 * t359 + (t447 + t578) * t355;
t374 = -qJD(1) * t45 + qJD(2) * t38 - qJD(3) * t213;
t317 = t335 * qJ(4);
t372 = t431 / 0.2e1 + t473 / 0.2e1;
t76 = t372 - t392;
t373 = qJD(1) * t76 + qJD(2) * t123 - qJD(3) * t317;
t369 = t422 * qJD(3) + t497;
t334 = qJD(3) * t441;
t331 = t351 * pkin(7) * t525;
t329 = t349 * t457;
t328 = t349 * t341;
t278 = t293 * t531;
t175 = (t399 + t436) * t353;
t172 = (t400 + t432) * t353;
t164 = t298 * t565 + t349 * t434 + t355 * t436 + t569;
t161 = t349 * t437 + t295 * t565 + t355 * t432 + t516 / 0.2e1;
t124 = t477 + t389;
t96 = t97 * qJD(5);
t81 = t226 * t438 + t227 * t445 + t424;
t78 = t376 + t386;
t77 = t372 + t392;
t60 = t227 * t562 + t294 * t438 + t353 * t424;
t59 = t226 * t562 + t294 * t444 + t355 * t425;
t51 = t360 * t391 + t390;
t48 = t141 * t565 + t226 * t443 + t355 * t447 + t535 / 0.2e1;
t47 = t140 * t566 + t356 * t452 + t355 * t446 - t536 / 0.2e1;
t44 = t141 * t445 + t454 - t544 / 0.2e1 + t359 * t424;
t42 = t140 * t444 + t456 - t543 / 0.2e1 + t356 * t425;
t39 = t364 + t388;
t37 = t365 + t387;
t28 = t366 - t394;
t26 = t368 - t396;
t24 = t453 + t549 / 0.2e1 + t267 * t443 + t370 * t353;
t23 = t455 - t550 / 0.2e1 + t356 * t449 + t371 * t353;
t14 = t542 * t575 + t269 * t541 / 0.2e1 + t293 * t477 + (t557 / 0.2e1 + pkin(3) * t440) * t357 + t378 - t393;
t7 = t362 - t397;
t6 = t379 + t367;
t1 = qJ(4) * t449 + t353 * t448 + t363 - t395;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t50 + qJD(3) * t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t31 + qJD(3) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t354 * t501, -t458, 0, 0, 0, 0, 0, 0, 0, 0, (-t358 * t500 - t361 * t482) * t354, (t357 * t501 - t361 * t481) * t354, (t351 + t352) * t458, t495 + (t331 + (pkin(7) * t352 * t361 - pkin(2) * t358) * t354) * qJD(2), 0, 0, 0, 0, 0, 0, (t267 * t360 + t351 * t472) * qJD(2) + t59 * qJD(3), (t268 * t360 + t351 * t471) * qJD(2) + t60 * qJD(3), t51 * qJD(3) + (t267 * t355 - t268 * t353) * t502, t510 + (-t267 * t269 + t268 * t270 + t331) * qJD(2) + t14 * qJD(3) + t81 * qJD(4), 0, 0, 0, 0, 0, 0, (t196 * t528 + t267 * t295) * qJD(2) + t23 * qJD(3) + t44 * qJD(5), (-t197 * t528 + t267 * t298) * qJD(2) + t24 * qJD(3) + t42 * qJD(5), (-t196 * t298 - t197 * t295) * qJD(2) + t6 * qJD(3), t532 + (-t142 * t196 + t143 * t197 + t241 * t267) * qJD(2) + t1 * qJD(3) + t28 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t357 * t458 - t489, -t360 * t458 + t490, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t59 - t294 * t483, qJD(2) * t60 + t353 * t489, t51 * qJD(2) - t335 * t490, t511 + t14 * qJD(2) + (-qJ(4) * t293 * t350 - pkin(3) * t294 - t278) * qJD(3) + t77 * qJD(4), 0, 0, 0, 0, 0, 0, t23 * qJD(2) + (-t184 * t355 - t293 * t530) * qJD(3) + t48 * qJD(5), t24 * qJD(2) + (t185 * t355 - t293 * t529) * qJD(3) + t47 * qJD(5), t6 * qJD(2) + (-t551 - t552) * t499, t533 + t1 * qJD(2) + (-t184 * t265 + t185 * t266 - t278) * qJD(3) + t26 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t81 + qJD(3) * t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t28 + qJD(3) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t44 + qJD(3) * t48 - qJD(5) * t141, qJD(2) * t42 + qJD(3) * t47 + qJD(5) * t140, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t495, 0, 0, 0, 0, 0, 0, -t58 * qJD(3), -t57 * qJD(3), t52 * qJD(3), -qJD(3) * t13 - qJD(4) * t80 - t510, 0, 0, 0, 0, 0, 0, qJD(3) * t22 + qJD(5) * t43, qJD(3) * t25 + qJD(5) * t41, -qJD(3) * t5, qJD(3) * t2 + qJD(4) * t29 - t532; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t457, t338 * qJD(3), 0, -t457, 0, 0, -pkin(2) * t482, -pkin(2) * t481, 0, 0, t350 * t457, -0.2e1 * t357 * t428, -t306 * qJD(3), t329, t305 * qJD(3), -t457, -qJD(3) * t118 + t355 * t462, qJD(3) * t119 - t353 * t462, -qJD(3) * t63 + qJD(4) * t302, qJD(3) * t74 - qJD(4) * t134, (qJD(3) * t299 - qJD(5) * t295) * t298, qJD(3) * t99 + qJD(5) * t151, qJD(3) * t170 - t295 * t461, (qJD(3) * t297 + qJD(5) * t298) * t295, qJD(3) * t169 - t298 * t461, t329, qJD(3) * t34 + qJD(4) * t219 + qJD(5) * t71, -qJD(3) * t35 + qJD(4) * t220 - qJD(5) * t70, qJD(3) * t20 - qJD(4) * t150, qJD(3) * t21 - qJD(4) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t341, t484, t481, -t341, -t482, 0, -pkin(7) * t481 - t479, pkin(7) * t482 - t478, 0, 0, (t350 * t502 + t466) * t360, (-0.2e1 * t353 * t467 + (-t349 + t350) * qJD(3)) * t360, t353 * t482 - t486, t328 - t428, t355 * t482 + t487, -t341, -t343 * qJD(3) + t353 * t369 - t410, t342 * qJD(3) + t355 * t369 - t409, qJD(3) * t411 - t412, (-pkin(3) * t556 + qJ(4) * t411) * qJD(3) + t124 * qJD(4) + t418, t581 + (t464 + t503) * t299, t508 - t96 + (-t297 * t359 - t299 * t356) * t499, t491 + (-t299 * t355 + t349 * t512) * qJD(3) + t161 * qJD(5), -t581 + (t465 + t504) * t297, t492 + (t297 * t355 - t349 * t515) * qJD(3) + t164 * qJD(5), t328 + (-t483 + qJD(5) / 0.2e1) * t527, -t146 * t483 + t172 * qJD(4) + t39 * qJD(5) + (qJ(4) * t297 + t243 * t356 - t265 * t360) * t499 + t417, t147 * t483 + t175 * qJD(4) + t37 * qJD(5) + (qJ(4) * t299 + t243 * t359 - t266 * t360) * t499 + t416, (t265 * t299 - t266 * t297 + (-t553 - t554) * t353) * qJD(3) + t78 * qJD(4) + t419, (qJ(4) * t243 * t353 - t146 * t265 + t147 * t266) * qJD(3) + t7 * qJD(4) + t420; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t467 + t499) * t360, -t383 * t360, t488, qJD(3) * t124 - t408, 0, 0, 0, 0, 0, 0, qJD(3) * t172 + t506, qJD(3) * t175 + t505, qJD(3) * t78 - t507, qJD(3) * t7 + t415; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t579, t407, t161 * qJD(3) + t295 * t421, t579, t164 * qJD(3) + t298 * t421, t334, qJD(3) * t39 - qJD(5) * t143 + t413, qJD(3) * t37 + qJD(5) * t142 + t414, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * qJD(2), t57 * qJD(2), -t52 * qJD(2), qJD(2) * t13 - qJD(4) * t76 - t511, 0, 0, 0, 0, 0, 0, -qJD(2) * t22 - qJD(5) * t45, -qJD(2) * t25 - qJD(5) * t46, qJD(2) * t5, -qJD(2) * t2 + qJD(4) * t27 - t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t341, -t484, 0, t341, 0, 0, t479, t478, 0, 0, -t350 * t341, 0.2e1 * t353 * t426, t486, -t328, -t487, t341, t410, t409, t412, -qJD(4) * t123 - t418, -t299 * t503 + t581, -t96 - t508, qJD(5) * t162 - t491, -t297 * t504 - t581, -qJD(5) * t163 - t492, qJD(5) * t442 - t328, qJD(4) * t173 + qJD(5) * t38 - t417, qJD(4) * t174 + qJD(5) * t36 - t416, qJD(4) * t79 - t419, qJD(4) * t8 - t420; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t335 * qJD(4), t317 * qJD(4), -t427, t303 * qJD(5), t353 * t460, t427, t353 * t459, 0, qJD(4) * t304 - qJD(5) * t213, qJD(4) * t307 + qJD(5) * t212, 0, t102 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t485, -t373, 0, 0, 0, 0, 0, 0, t402, t401, t509, -t380; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t580, -t406, t356 * t403 + t494, -t580, t359 * t403 - t493, qJD(2) * t442, -qJD(5) * t266 + t374, qJD(5) * t265 + t375, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t80 + qJD(3) * t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t29 - qJD(3) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, t353 * t341, -t488, qJD(3) * t123 + t408, 0, 0, 0, 0, 0, 0, -qJD(3) * t173 - t356 * t461 - t506, -qJD(3) * t174 - t359 * t461 - t505, -qJD(3) * t79 + t507, -qJD(3) * t8 - t415; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485, t373, 0, 0, 0, 0, 0, 0, -t402 + t460, -t401 + t459, -t509, t380; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t381 * t356, t381 * t359, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t43 + qJD(3) * t45, -qJD(2) * t41 + qJD(3) * t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t579, -t407, -qJD(3) * t162 + t295 * t468, -t579, qJD(3) * t163 + t298 * t468, t334, -qJD(3) * t38 + t356 * t463 - t413, -qJD(3) * t36 + t359 * t463 - t414, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t580, t406, -t355 * t465 - t494, t580, -t355 * t464 + t493, qJD(2) * t441, -t356 * t498 - t374, -t359 * t498 - t375, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t383 * t356, t383 * t359, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
