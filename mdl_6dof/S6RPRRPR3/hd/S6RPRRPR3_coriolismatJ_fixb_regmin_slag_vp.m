% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x29]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:18
% EndTime: 2019-03-09 05:07:36
% DurationCPUTime: 9.43s
% Computational Cost: add. (4835->552), mult. (9906->723), div. (0->0), fcn. (9328->8), ass. (0->434)
t469 = qJD(4) - qJD(6);
t347 = sin(qJ(6));
t348 = sin(qJ(4));
t349 = sin(qJ(3));
t350 = cos(qJ(6));
t535 = t350 * t349;
t351 = cos(qJ(4));
t541 = t349 * t351;
t222 = t347 * t541 - t348 * t535;
t517 = t222 * t469;
t271 = t347 * t348 + t350 * t351;
t224 = t271 * t349;
t352 = cos(qJ(3));
t473 = t352 * qJD(1);
t455 = t224 * t473;
t515 = t469 * t224;
t620 = t515 - t455;
t441 = t222 * t473;
t619 = t517 - t441;
t536 = t350 * t348;
t275 = -t347 * t351 + t536;
t508 = qJD(3) * t275;
t601 = t275 / 0.2e1;
t602 = t271 / 0.2e1;
t57 = t222 * t601 + t224 * t602;
t618 = t57 * qJD(1) + t271 * t508;
t509 = qJD(1) * t224;
t617 = qJD(3) * t57 + t222 * t509;
t543 = t348 * t349;
t324 = pkin(4) * t543;
t325 = sin(pkin(10)) * pkin(1) + pkin(7);
t534 = t351 * qJ(5);
t422 = -t325 + t534;
t164 = -t324 + (-t348 * pkin(5) + t422) * t349;
t613 = pkin(8) - pkin(9);
t298 = t613 * t348;
t299 = t613 * t351;
t421 = t347 * t298 + t350 * t299;
t593 = -t352 / 0.2e1;
t545 = t348 * qJ(5);
t609 = pkin(4) + pkin(5);
t610 = -t351 * t609 - t545;
t253 = pkin(3) - t610;
t604 = t253 / 0.2e1;
t364 = t164 * t601 + t224 * t604 + t421 * t593;
t616 = t469 * t421;
t41 = t271 * t222 - t224 * t275;
t60 = t222 ^ 2 - t224 ^ 2;
t402 = qJD(1) * t60 + qJD(3) * t41;
t110 = t271 ^ 2 - t275 ^ 2;
t397 = qJD(1) * t41 + qJD(3) * t110;
t113 = t350 * t298 - t347 * t299;
t615 = t469 * t113;
t614 = 0.2e1 * t348;
t612 = t271 * t469;
t416 = t469 * t275;
t528 = t352 * t325;
t281 = t351 * t528;
t326 = -cos(pkin(10)) * pkin(1) - pkin(2);
t589 = t349 * pkin(8);
t409 = -t352 * pkin(3) - t589;
t368 = t409 + t326;
t167 = t348 * t368 + t281;
t529 = t352 * qJ(5);
t144 = t167 - t529;
t608 = -t144 / 0.2e1;
t228 = t351 * t368;
t462 = t348 * t528;
t166 = -t228 + t462;
t607 = -t166 / 0.2e1;
t194 = -t349 * t422 + t324;
t606 = t194 / 0.2e1;
t588 = t351 * pkin(4);
t407 = t545 + t588;
t252 = t407 * t349;
t605 = -t252 / 0.2e1;
t263 = -t348 * t609 + t534;
t603 = t263 / 0.2e1;
t280 = t325 * t543;
t600 = -t280 / 0.2e1;
t591 = pkin(4) * t348;
t297 = -t534 + t591;
t599 = t297 / 0.2e1;
t587 = t352 * pkin(8);
t590 = t349 * pkin(3);
t300 = -t587 + t590;
t598 = -t300 / 0.2e1;
t597 = t348 / 0.2e1;
t596 = -t349 / 0.2e1;
t595 = t349 / 0.2e1;
t594 = t351 / 0.2e1;
t592 = t352 / 0.2e1;
t341 = t352 * pkin(4);
t586 = t469 * t41;
t585 = t469 * t57;
t165 = (-t325 + t263) * t352;
t304 = t352 * t536;
t530 = t351 * t352;
t225 = t347 * t530 - t304;
t323 = pkin(9) * t543;
t99 = t144 + t323;
t582 = t347 * t99;
t467 = pkin(9) * t541;
t551 = t325 * t348;
t358 = -t467 - t228 + t341 + (pkin(5) + t551) * t352;
t69 = t350 * t358;
t42 = -t69 + t582;
t340 = t349 * pkin(4);
t106 = -t349 * pkin(5) - t280 - t340 + (-pkin(9) * t352 - t300) * t351;
t540 = t350 * t106;
t277 = t348 * t300;
t338 = t349 * qJ(5);
t175 = -t325 * t541 + t277 + t338;
t542 = t348 * t352;
t126 = pkin(9) * t542 + t175;
t546 = t347 * t126;
t2 = (t540 - t546) * t352 + t42 * t349 + t165 * t222 + t164 * t225;
t584 = t2 * qJD(1);
t226 = t352 * t271;
t357 = t347 * t358;
t581 = t350 * t99;
t43 = t357 + t581;
t537 = t350 * t126;
t549 = t347 * t106;
t3 = (t537 + t549) * t352 - t43 * t349 - t165 * t224 - t164 * t226;
t583 = t3 * qJD(1);
t561 = t252 * t348;
t566 = t194 * t351;
t570 = t167 * t352;
t45 = t570 + (t561 + t566) * t349;
t580 = qJD(1) * t45;
t344 = t349 ^ 2;
t571 = t166 * t352;
t74 = -t344 * t551 - t571;
t579 = qJD(1) * t74;
t531 = t351 * t344;
t75 = -t325 * t531 - t570;
t578 = qJD(1) * t75;
t283 = t350 * qJ(5) - t347 * t609;
t354 = t283 * t592 - t357 / 0.2e1;
t122 = -t166 + t467;
t548 = t347 * t122;
t438 = -t548 / 0.2e1;
t123 = t323 + t167;
t461 = t123 / 0.2e1 - t99 / 0.2e1;
t11 = t350 * t461 + t354 + t438;
t577 = t11 * qJD(1);
t282 = t347 * qJ(5) + t350 * t609;
t410 = -t69 / 0.2e1 + t282 * t593;
t539 = t350 * t122;
t432 = -t539 / 0.2e1;
t12 = -t347 * t461 + t410 + t432;
t576 = t12 * qJD(1);
t575 = t144 * t351;
t574 = t144 * t352;
t573 = t164 * t222;
t572 = t164 * t224;
t569 = t175 * t351;
t532 = t351 * t300;
t420 = -t280 - t532;
t178 = -t340 + t420;
t568 = t178 * t348;
t567 = t194 * t348;
t145 = t166 + t341;
t425 = t607 + t145 / 0.2e1;
t426 = t608 + t167 / 0.2e1;
t359 = t348 * t426 + t351 * t425;
t21 = t252 * t593 + t349 * t359;
t565 = t21 * qJD(1);
t205 = t610 * t349;
t538 = t350 * t123;
t49 = -t538 + t548;
t22 = t205 * t222 - t49 * t352 - t572;
t564 = t22 * qJD(1);
t547 = t347 * t123;
t50 = t539 + t547;
t23 = -t205 * t224 + t50 * t352 - t573;
t563 = t23 * qJD(1);
t25 = (-t341 / 0.2e1 - t425) * t351 + (-t529 / 0.2e1 - t426) * t348;
t562 = t25 * qJD(1);
t560 = t253 * t271;
t27 = -t167 * t541 + (t575 + (t145 - t166) * t348) * t349;
t559 = t27 * qJD(1);
t558 = t275 * t352;
t28 = -t145 * t530 - t178 * t541 + (t175 * t349 + t574) * t348;
t557 = t28 * qJD(1);
t284 = -pkin(3) - t407;
t556 = t284 * t348;
t29 = -t42 * t352 + t573;
t555 = t29 * qJD(1);
t30 = -t43 * t352 + t572;
t554 = t30 * qJD(1);
t195 = (t297 + t325) * t352;
t400 = t194 * t352 + t195 * t349;
t31 = -t144 * t349 + t175 * t352 + t351 * t400;
t553 = t31 * qJD(1);
t32 = -t145 * t349 + t178 * t352 + t348 * t400;
t552 = t32 * qJD(1);
t343 = t348 ^ 2;
t550 = t343 * t352;
t544 = t348 * t224;
t533 = t351 * t275;
t46 = -t194 * t543 + t252 * t541 - t571;
t527 = t46 * qJD(1);
t51 = t166 * t349 + (-t420 - 0.2e1 * t280) * t352;
t526 = t51 * qJD(1);
t52 = t277 * t352 + (-t167 + t281) * t349;
t525 = t52 * qJD(1);
t53 = -t222 * t226 - t224 * t225;
t524 = t53 * qJD(1);
t54 = t194 * t541 + t574;
t523 = t54 * qJD(1);
t89 = t349 * t222 - t352 * t225;
t522 = t89 * qJD(1);
t90 = -t224 * t349 + t226 * t352;
t521 = t90 * qJD(1);
t427 = t226 / 0.2e1;
t429 = t530 / 0.2e1;
t435 = t542 / 0.2e1;
t514 = t347 * t435 + t350 * t429;
t137 = t427 + t514;
t428 = -t226 / 0.2e1;
t430 = -t530 / 0.2e1;
t436 = -t542 / 0.2e1;
t513 = t347 * t436 + t350 * t430;
t138 = t428 + t513;
t520 = -t138 * qJD(4) - t137 * qJD(6);
t136 = t428 + t514;
t139 = t427 + t513;
t519 = -t136 * qJD(4) - t139 * qJD(6);
t512 = t304 / 0.2e1 + t347 * t430;
t511 = -t304 / 0.2e1 + t347 * t429;
t345 = t351 ^ 2;
t510 = t343 + t345;
t312 = t345 - t343;
t346 = t352 ^ 2;
t313 = t346 - t344;
t507 = qJD(3) * t348;
t506 = qJD(3) * t351;
t505 = qJD(4) * t348;
t504 = qJD(4) * t351;
t503 = qJD(6) * t253;
t440 = -t558 / 0.2e1;
t132 = t440 + t511;
t502 = t132 * qJD(1);
t501 = t132 * qJD(3);
t439 = t558 / 0.2e1;
t133 = t439 + t511;
t500 = t133 * qJD(3);
t134 = t440 + t512;
t499 = t134 * qJD(3);
t135 = t439 + t512;
t498 = t135 * qJD(1);
t497 = t135 * qJD(3);
t496 = t136 * qJD(3);
t495 = t137 * qJD(1);
t494 = t137 * qJD(3);
t493 = t138 * qJD(1);
t492 = t138 * qJD(3);
t491 = t139 * qJD(3);
t490 = t166 * qJD(4);
t168 = t222 * t541 + t347 * t346;
t489 = t168 * qJD(1);
t169 = t224 * t541 + t350 * t346;
t488 = t169 * qJD(1);
t274 = t313 * t348;
t483 = t274 * qJD(1);
t276 = t351 * t346 - t531;
t482 = t276 * qJD(1);
t481 = t313 * qJD(1);
t480 = t347 * qJD(5);
t479 = t348 * qJD(5);
t478 = t349 * qJD(1);
t477 = t349 * qJD(3);
t476 = t349 * qJD(4);
t475 = t350 * qJD(5);
t474 = t351 * qJD(5);
t472 = t352 * qJD(3);
t471 = t352 * qJD(4);
t470 = t352 * qJD(5);
t468 = t600 - t340;
t466 = pkin(8) * t505;
t465 = pkin(8) * t504;
t464 = pkin(3) * t597;
t463 = t587 / 0.2e1;
t459 = t351 * t478;
t458 = t271 * t477;
t457 = t348 * t471;
t456 = t351 * t471;
t454 = t275 * t477;
t452 = t326 * t478;
t451 = t326 * t473;
t450 = t347 * t470;
t449 = t348 * t504;
t318 = t348 * t506;
t448 = t348 * t474;
t447 = t349 * t472;
t446 = t349 * t479;
t445 = t349 * t473;
t444 = t350 * t470;
t443 = t351 * t477;
t442 = t349 * t474;
t437 = t222 * t597;
t434 = t284 * t595;
t433 = -t541 / 0.2e1;
t431 = t271 * t594;
t424 = t600 - t340 / 0.2e1;
t423 = t343 / 0.2e1 - t345 / 0.2e1;
t233 = (-0.1e1 / 0.2e1 + t423) * t349;
t419 = t233 * qJD(1) - t318;
t255 = t423 * t349;
t418 = t255 * qJD(1) - t318;
t295 = t348 * qJD(1) * t531;
t417 = t255 * qJD(3) + t295;
t415 = t469 * t352;
t322 = -qJD(4) + t473;
t414 = qJD(6) + t473;
t413 = t348 * t459;
t412 = t348 * t443;
t411 = t473 - qJD(4) / 0.2e1;
t408 = t297 * t595 + t606;
t406 = -t284 * t352 + t589;
t10 = (t575 / 0.2e1 - t195 / 0.2e1 + t145 * t597) * t352 + (t569 / 0.2e1 + t606 + t568 / 0.2e1) * t349;
t20 = t144 * t175 + t145 * t178 + t194 * t195;
t405 = t20 * qJD(1) + t10 * qJD(2);
t24 = -t144 * t166 + t145 * t167 + t194 * t252;
t404 = t24 * qJD(1) + t21 * qJD(2);
t401 = t568 + t569;
t196 = (-0.1e1 + t510) * t352 * t349;
t399 = t10 * qJD(1) + t196 * qJD(2);
t185 = t284 * t351 + t297 * t348;
t259 = -t277 / 0.2e1;
t306 = pkin(8) * t435;
t34 = t306 - t561 / 0.2e1 - t566 / 0.2e1 - t338 + t259 + (t556 / 0.2e1 + (-t297 / 0.2e1 + t325 / 0.2e1) * t351) * t349;
t396 = -t34 * qJD(1) + t185 * qJD(3);
t186 = t297 * t351 - t556;
t390 = t434 + t463;
t372 = t605 + t390;
t387 = t408 * t348;
t36 = t387 + (t598 + t372) * t351 + t468;
t395 = -t36 * qJD(1) + t186 * qJD(3);
t394 = t322 * t349;
t393 = t463 - t590 / 0.2e1;
t392 = -t591 / 0.2e1 + t534 / 0.2e1;
t391 = -t175 * qJ(5) / 0.2e1 + t178 * pkin(4) / 0.2e1;
t258 = t277 / 0.2e1;
t307 = pkin(8) * t436;
t179 = t349 * t464 + t258 + t307;
t389 = pkin(3) * t506 - t179 * qJD(1);
t180 = (t598 + t393) * t351;
t388 = pkin(3) * t507 - t180 * qJD(1);
t386 = t164 * t602 + t222 * t604;
t384 = -t549 / 0.2e1 - t537 / 0.2e1;
t383 = -t546 / 0.2e1 + t540 / 0.2e1;
t71 = t437 + (t431 + t350 / 0.2e1) * t349;
t378 = t71 * qJD(1) + t271 * t507;
t70 = -t544 / 0.2e1 + (t347 / 0.2e1 - t533 / 0.2e1) * t349;
t377 = -t70 * qJD(1) + t275 * t507;
t47 = t567 / 0.2e1 + (t598 + t390) * t351 + t424;
t376 = t47 * qJD(1) + t284 * t507;
t375 = t351 * t394;
t273 = t312 * t344;
t374 = t273 * qJD(1) + 0.2e1 * t412;
t373 = -t312 * qJD(3) + 0.2e1 * t413;
t355 = t205 * t602 + t222 * t603 - t364;
t363 = t282 * t596 + t383;
t4 = t355 + t363;
t61 = -t253 * t275 + t263 * t271;
t371 = -qJD(1) * t4 + qJD(2) * t133 - qJD(3) * t61;
t62 = t263 * t275 + t560;
t356 = -t113 * t593 + t205 * t601 + t224 * t603 + t386;
t362 = t283 * t596 + t384;
t7 = t356 + t362;
t370 = -qJD(1) * t7 + qJD(2) * t136 - qJD(3) * t62;
t176 = (t599 + t392) * t352;
t353 = t359 * pkin(8) + t194 * t599 + t252 * t284 / 0.2e1;
t9 = t353 + t391;
t369 = t284 * t297 * qJD(3) + t9 * qJD(1) - t176 * qJD(2);
t15 = -t364 + t383;
t367 = t15 * qJD(1) + t134 * qJD(2) - t253 * t508;
t365 = -t113 * t592 - t386;
t16 = -t365 + t384;
t366 = t16 * qJD(1) + t139 * qJD(2) + qJD(3) * t560;
t361 = -qJD(4) * t407 + t474;
t286 = t345 * t344 + t346;
t360 = t286 * qJD(1) + t412 - t471;
t337 = t350 * qJD(4);
t336 = t347 * qJD(4);
t332 = t345 * t352;
t331 = -t478 / 0.2e1;
t330 = t478 / 0.2e1;
t329 = -t477 / 0.2e1;
t328 = t477 / 0.2e1;
t321 = t351 * t472;
t320 = t351 * t473;
t319 = t348 * t477;
t296 = t348 * t442;
t285 = t322 * qJ(5);
t270 = -t320 + t504;
t269 = t350 * t473 - t337;
t268 = t347 * t473 - t336;
t265 = t411 * t349;
t262 = t532 / 0.2e1;
t250 = -t319 + t456;
t249 = t348 * t472 + t351 * t476;
t248 = -t443 - t457;
t247 = t348 * t476 - t321;
t246 = t343 * qJD(3) + t413;
t241 = t255 * qJD(4);
t234 = t343 * t596 + t345 * t595 + t596;
t232 = t350 * t414 - t337;
t231 = t347 * t414 - t336;
t221 = -t348 * t445 + t321;
t220 = t348 * t394;
t219 = (qJD(6) / 0.2e1 + t411) * t349;
t177 = t297 * t593 + t352 * t392;
t154 = t167 * qJD(4);
t103 = t351 * t393 + t262 + t280;
t102 = t307 + t259 + (t325 * t351 + t464) * t349;
t73 = t544 / 0.2e1 + (t533 + t347) * t595;
t72 = t437 + t349 * t431 - t535 / 0.2e1;
t48 = pkin(8) * t430 + t284 * t433 - t567 / 0.2e1 - t532 / 0.2e1 + t424;
t35 = t351 * t372 + t262 + t387 - t468;
t33 = t306 + t338 + t325 * t433 + t258 - t408 * t351 + (t434 + t605) * t348;
t26 = t351 * t607 + t348 * t608 + t167 * t597 + t145 * t594 + (-t545 / 0.2e1 - t588 / 0.2e1) * t352;
t18 = t364 + t383;
t17 = t365 + t384;
t14 = t581 / 0.2e1 + t438 + t538 / 0.2e1 - t354;
t13 = -t582 / 0.2e1 + t432 - t547 / 0.2e1 - t410;
t8 = t353 - t391;
t6 = t356 - t362;
t5 = t355 - t363;
t1 = qJD(3) * t10 + qJD(4) * t21;
t19 = [0, 0, 0, 0, t447, t313 * qJD(3), 0, 0, 0, t326 * t477, t326 * t472, -t344 * t449 + t345 * t447, -t273 * qJD(4) - 0.2e1 * t352 * t412, -t276 * qJD(3) + t349 * t457, t274 * qJD(3) + t349 * t456, -t447, -qJD(3) * t51 - qJD(4) * t75, qJD(3) * t52 + qJD(4) * t74, t32 * qJD(3) + t45 * qJD(4) - t344 * t448, -t28 * qJD(3) - t27 * qJD(4) + t352 * t446, -t31 * qJD(3) - t46 * qJD(4) + t286 * qJD(5), qJD(3) * t20 + qJD(4) * t24 - qJD(5) * t54 (qJD(3) * t226 + t517) * t224, qJD(3) * t53 - t469 * t60, t90 * qJD(3) + t222 * t415, t89 * qJD(3) + t224 * t415, -t447, qJD(3) * t2 + qJD(4) * t22 + qJD(5) * t168 + qJD(6) * t30, -qJD(3) * t3 - qJD(4) * t23 + qJD(5) * t169 - qJD(6) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t445, t481, t472, -t477, 0, -t325 * t472 + t452, t325 * t477 + t451, -t241 + (t345 * t478 + t318) * t352 (t332 - t550) * qJD(3) + (-qJD(4) - t473) * t541 * t614, t319 - t482, t443 + t483, -t265, -t526 + (t348 * t409 - t281) * qJD(3) + t103 * qJD(4), t525 + (t351 * t409 + t462) * qJD(3) + t102 * qJD(4), t552 + (-t195 * t351 - t348 * t406) * qJD(3) + t35 * qJD(4) + t234 * qJD(5), qJD(3) * t401 + t26 * qJD(4) - t557, -t553 + (-t195 * t348 + t351 * t406) * qJD(3) + t33 * qJD(4) + t296 (pkin(8) * t401 + t195 * t284) * qJD(3) + t8 * qJD(4) + t48 * qJD(5) + t405 (t508 + t509) * t226 + t585, t524 + (-t275 * t225 - t226 * t271) * qJD(3) - t586, -t454 + t519 + t521, -t134 * qJD(4) - t133 * qJD(6) + t458 + t522, -t219, t584 + (-t113 * t349 + t165 * t271 + t253 * t225) * qJD(3) + t5 * qJD(4) + t72 * qJD(5) + t18 * qJD(6), -t583 + (t165 * t275 + t253 * t226 + t349 * t421) * qJD(3) + t6 * qJD(4) + t73 * qJD(5) + t17 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t417, -t374, t220, t375, t328, qJD(3) * t103 - t154 - t578, qJD(3) * t102 + t490 + t579, qJD(3) * t35 - t154 + t580, -t559 + t26 * qJD(3) + (-t349 * t534 + t324) * qJD(4) - t446, t33 * qJD(3) - t470 - t490 - t527, t8 * qJD(3) + (-pkin(4) * t167 - qJ(5) * t166) * qJD(4) + t144 * qJD(5) + t404, t617, -t402, -t496 - t619, -t499 - t620, t328, t5 * qJD(3) + t49 * qJD(4) + t14 * qJD(6) - t450 + t564, t6 * qJD(3) + t50 * qJD(4) + t13 * qJD(6) - t444 - t563; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234 * qJD(3) - t295, t220, t360, qJD(3) * t48 + qJD(4) * t144 - t523, 0, 0, 0, 0, 0, t72 * qJD(3) - t347 * t471 + t489, t73 * qJD(3) - t350 * t471 + t488; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t617, t402, -t491 + t619, -t500 + t620, t329, qJD(3) * t18 + qJD(4) * t14 - qJD(6) * t43 + t554, qJD(3) * t17 + qJD(4) * t13 + qJD(6) * t42 - t555; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t477, -t472, 0, 0, 0, 0, 0, t248, -t250, t248 (t332 + t550) * qJD(3), t250, t284 * t477 + t177 * qJD(4) + (pkin(8) * qJD(3) * t510 + t479) * t352 + t399, 0, 0, 0, 0, 0, -t135 * qJD(4) - t132 * qJD(6) - t458, -t454 + t520; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t249, t247, -t249, 0, -t247, t177 * qJD(3) - t252 * qJD(4) + t442 + t565, 0, 0, 0, 0, 0, -t497 - t515, -t492 + t517; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t501 + t515, -t494 - t517; 0, 0, 0, 0, -t445, -t481, 0, 0, 0, -t452, -t451, -t345 * t445 - t241, t375 * t614, -t456 + t482, t457 - t483, t265, qJD(4) * t180 + t526, qJD(4) * t179 - t525, qJD(4) * t36 - qJD(5) * t233 - t552, -t25 * qJD(4) - t351 * t470 + t557, t34 * qJD(4) + t296 + t553, qJD(4) * t9 - qJD(5) * t47 - t405, -t226 * t509 + t585, -t524 - t586, t520 - t521, -qJD(4) * t132 - qJD(6) * t135 - t522, t219, qJD(4) * t4 + qJD(5) * t71 - qJD(6) * t15 - t584, qJD(4) * t7 - qJD(5) * t70 - qJD(6) * t16 + t583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t176 - t399, 0, 0, 0, 0, 0, -qJD(4) * t133 - qJD(6) * t134, t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t449, t312 * qJD(4), 0, 0, 0, -pkin(3) * t505, -pkin(3) * t504, -t186 * qJD(4) + t448, 0, -t185 * qJD(4) + t343 * qJD(5) (qJD(4) * t297 - t479) * t284, t271 * t416, -t469 * t110, 0, 0, 0, t61 * qJD(4) + t271 * t479 + t275 * t503, t62 * qJD(4) - t271 * t503 + t275 * t479; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t418, -t373, t270, t322 * t348, t331, -t388 - t465, -t389 + t466, -t395 - t465, t361 - t562, -t396 - t466, pkin(8) * t361 + t369, t618, -t397, -t493 - t612, -t502 - t416, t331, -t371 - t616, -t370 - t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t419, t270, t246, -t376 + t465, 0, 0, 0, 0, 0, t378, t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t618, t397, -t495 + t612, t416 - t498, t330, -t367 + t616, -t366 + t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t417, t374, t221 (-t459 - t507) * t352, t328, -qJD(3) * t180 + t578, -qJD(3) * t179 - t579, -qJD(3) * t36 - t580, qJD(3) * t25 + t559, -t34 * qJD(3) - t470 + t527, -qJ(5) * t470 - t9 * qJD(3) - t404, -t617, t402, -t441 + t492, -t455 + t501, t328, -t4 * qJD(3) - t11 * qJD(6) - t450 - t564, -t7 * qJD(3) - t12 * qJD(6) - t444 + t563; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t176 - t565, 0, 0, 0, 0, 0, t500, t496; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t418, t373, t320, -t348 * t473, t330, t388, t389, t395, t562, t396, -t369, -t618, t397, t493, t502, t330, t371, t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5), 0, 0, 0, 0, 0, t283 * qJD(6) + t480, -t282 * qJD(6) + t475; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, -t285, 0, 0, 0, 0, 0, -t268, -t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283 * t469 - t577, -t282 * t469 - t576; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233 * qJD(3) + t295, t221, -t360, qJ(5) * t471 + t47 * qJD(3) + t523, 0, 0, 0, 0, 0, -t71 * qJD(3) + t347 * t415 - t489, t70 * qJD(3) + t350 * t415 - t488; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419, t320, -t246, t376, 0, 0, 0, 0, 0, -t378, -t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, t285, 0, 0, 0, 0, 0, t231, t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, -t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t617, -t402, t441 + t494, t455 + t497, t329, t15 * qJD(3) + t11 * qJD(4) + t450 - t554, t16 * qJD(3) + t12 * qJD(4) + t444 + t555; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t499, t491; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t618, -t397, t495, t498, t331, t367, t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t283 * qJD(4) - t480 + t577, t282 * qJD(4) - t475 + t576; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t19;
