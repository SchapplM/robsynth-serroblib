% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPP7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:15
% EndTime: 2019-03-09 04:52:30
% DurationCPUTime: 8.34s
% Computational Cost: add. (4824->632), mult. (8895->744), div. (0->0), fcn. (7318->4), ass. (0->436)
t335 = sin(qJ(4));
t323 = t335 * qJ(5);
t337 = cos(qJ(4));
t556 = t337 * pkin(4);
t372 = t323 + t556;
t354 = t372 * qJD(4);
t445 = t337 * qJD(5);
t591 = -t354 + t445;
t580 = pkin(4) + pkin(5);
t489 = t337 * t580;
t590 = -t489 - t323;
t553 = pkin(8) - qJ(6);
t336 = sin(qJ(3));
t557 = t336 * pkin(3);
t338 = cos(qJ(3));
t558 = pkin(8) * t338;
t380 = t557 - t558;
t363 = qJ(2) + t380;
t231 = t337 * t363;
t340 = -pkin(1) - pkin(7);
t493 = t336 * t340;
t161 = t335 * t493 - t231;
t490 = t337 * t338;
t432 = qJ(6) * t490;
t115 = -t161 + t432;
t496 = t335 * t340;
t393 = -pkin(4) + t496;
t73 = -t432 - t231 + (-pkin(5) + t393) * t336;
t589 = t115 + t73;
t488 = t337 * t340;
t433 = t336 * t488;
t162 = t335 * t363 + t433;
t324 = t336 * qJ(5);
t138 = t324 + t162;
t588 = (t138 / 0.2e1 - t162 / 0.2e1) * t335;
t332 = t335 ^ 2;
t333 = t337 ^ 2;
t281 = t333 + t332;
t443 = t338 * qJD(3);
t296 = t335 * t443;
t462 = qJD(4) * t337;
t302 = t336 * t462;
t583 = t336 ^ 2;
t441 = t583 / 0.2e1;
t334 = t338 ^ 2;
t313 = t334 * t337;
t571 = t313 / 0.2e1;
t202 = t571 + (t441 + 0.1e1 / 0.2e1) * t337;
t455 = t202 * qJD(1);
t86 = t455 + t296 + t302;
t300 = t337 * t443;
t463 = qJD(4) * t336;
t423 = t335 * t463;
t391 = t441 + t334 / 0.2e1;
t199 = (0.1e1 / 0.2e1 + t391) * t335;
t456 = t199 * qJD(1);
t587 = t456 - t300 + t423;
t485 = t338 * t340;
t267 = t337 * t485;
t325 = t338 * qJ(5);
t585 = t325 + t267 / 0.2e1;
t492 = t337 * qJ(5);
t230 = -t335 * t580 + t492;
t282 = t333 - t332;
t444 = t338 * qJD(1);
t424 = t337 * t444;
t385 = t335 * t424;
t351 = -qJD(3) * t282 + 0.2e1 * t385;
t582 = pkin(4) / 0.2e1;
t116 = t433 + (-t338 * t553 + qJ(2) + t557) * t335;
t92 = t324 + t116;
t581 = -t92 / 0.2e1;
t307 = t337 * t325;
t497 = t335 * t338;
t309 = pkin(4) * t497;
t474 = t307 - t309;
t137 = (-t335 * pkin(5) + t340) * t338 + t474;
t579 = -t137 / 0.2e1;
t166 = t590 * t338;
t578 = -t166 / 0.2e1;
t170 = -t474 - t485;
t577 = t170 / 0.2e1;
t217 = t372 * t338;
t576 = -t217 / 0.2e1;
t255 = t553 * t335;
t575 = t255 / 0.2e1;
t555 = t338 * pkin(3);
t258 = pkin(8) * t336 + t555;
t573 = -t258 / 0.2e1;
t266 = t335 * t485;
t572 = t266 / 0.2e1;
t570 = -t323 / 0.2e1;
t569 = t332 / 0.2e1;
t568 = -t335 / 0.2e1;
t567 = t335 / 0.2e1;
t566 = -t336 / 0.2e1;
t565 = t336 / 0.2e1;
t564 = -t337 / 0.2e1;
t563 = t337 / 0.2e1;
t562 = -t338 / 0.2e1;
t561 = t338 / 0.2e1;
t560 = t580 / 0.2e1;
t559 = -t580 / 0.2e1;
t554 = t338 * pkin(5);
t136 = (-t230 - t340) * t336;
t257 = t553 * t337;
t326 = t338 * pkin(4);
t495 = t336 * qJ(6);
t79 = -t554 + t266 - t326 + (-t258 + t495) * t337;
t549 = t79 * t336;
t242 = t335 * t258;
t148 = t267 + t242 + t325;
t96 = -t335 * t495 + t148;
t4 = t136 * t561 + t137 * t566 + (t561 * t92 + t565 * t96 + t575) * t337 + (t549 / 0.2e1 + t73 * t561 - t257 / 0.2e1) * t335;
t552 = t4 * qJD(1);
t7 = t136 * t137 + t73 * t79 + t92 * t96;
t551 = t7 * qJD(1);
t550 = t73 * t337;
t498 = t335 * t336;
t412 = -t498 / 0.2e1;
t261 = qJ(5) * t412;
t430 = t92 / 0.2e1 - t116 / 0.2e1;
t383 = t430 * t335;
t399 = t336 * t559;
t431 = t115 / 0.2e1 + t73 / 0.2e1;
t8 = t261 - t383 + (t399 + t431) * t337;
t548 = t8 * qJD(1);
t547 = t92 * t335;
t546 = t92 * t336;
t506 = t257 * t336;
t413 = -t506 / 0.2e1;
t529 = t137 * t335;
t357 = -t529 / 0.2e1 + t413;
t438 = -t326 / 0.2e1 + t572;
t346 = t357 + t438;
t408 = t495 / 0.2e1;
t220 = pkin(3) - t590;
t414 = t220 * t562;
t356 = t414 + t408;
t435 = -t554 / 0.2e1;
t39 = t435 + (t573 + t356) * t337 + t346;
t545 = qJD(1) * t39;
t41 = t115 * t336 + (t166 * t337 - t529) * t338;
t544 = qJD(1) * t41;
t518 = t166 * t335;
t528 = t137 * t337;
t42 = t116 * t336 + (t518 + t528) * t338;
t543 = qJD(1) * t42;
t44 = t547 - t550;
t43 = t44 * t338;
t542 = qJD(1) * t43;
t541 = qJD(1) * t44;
t515 = t170 * t337;
t519 = t162 * t336;
t47 = -t519 + (t217 * t335 + t515) * t338;
t540 = qJD(1) * t47;
t512 = t217 * t337;
t516 = t170 * t335;
t520 = t161 * t336;
t48 = -t520 + (-t512 + t516) * t338;
t539 = qJD(1) * t48;
t49 = t137 * t490 + t546;
t538 = qJD(1) * t49;
t142 = t336 * t393 - t231;
t524 = t142 * t337;
t527 = t138 * t335;
t55 = -t524 + t527;
t537 = qJD(1) * t55;
t56 = t138 * t336 - t170 * t490;
t536 = qJD(1) * t56;
t83 = -t334 * t496 - t520;
t535 = qJD(1) * t83;
t84 = -t334 * t488 - t519;
t534 = qJD(1) * t84;
t342 = (t337 * t431 - t383) * t336 + t166 * t561;
t11 = t570 - t489 / 0.2e1 + t342;
t533 = t11 * qJD(1);
t13 = t115 * t92 + t116 * t73 + t137 * t166;
t532 = t13 * qJD(1);
t531 = t136 * t335;
t530 = t136 * t337;
t526 = t138 * t337;
t256 = pkin(4) * t335 - t492;
t169 = (-t256 + t340) * t336;
t491 = t337 * t258;
t392 = t266 - t491;
t152 = -t326 + t392;
t521 = t152 * t335;
t523 = t148 * t337;
t14 = (t526 / 0.2e1 - t169 / 0.2e1 + t142 * t567) * t338 + (t523 / 0.2e1 + t577 + t521 / 0.2e1) * t336;
t525 = t14 * qJD(1);
t15 = (t73 * t336 - t338 * t79) * t337 + (t338 * t96 - t546) * t335;
t522 = t15 * qJD(1);
t517 = t169 * t337;
t18 = ((-t116 + t92) * t337 + t589 * t335) * t338;
t514 = t18 * qJD(1);
t396 = -t161 / 0.2e1 + t142 / 0.2e1;
t379 = t396 * t337;
t415 = t162 * t567;
t98 = t138 * t498;
t343 = (t415 + t379) * t336 - t98 / 0.2e1 + t217 * t562;
t20 = t570 - t556 / 0.2e1 + t343;
t513 = t20 * qJD(1);
t22 = t549 - t137 * t498 + (t73 + t531) * t338;
t511 = t22 * qJD(1);
t510 = t230 * t335;
t245 = -pkin(3) - t372;
t509 = t245 * t336;
t508 = t255 * t336;
t507 = t255 * t338;
t505 = t257 * t338;
t27 = (t92 + t530) * t338 + (t96 - t528) * t336;
t504 = t27 * qJD(1);
t28 = -t138 * t161 + t142 * t162 + t170 * t217;
t503 = t28 * qJD(1);
t409 = t324 / 0.2e1;
t260 = t335 * t409;
t29 = t260 + t588 + (pkin(4) * t565 - t396) * t337;
t502 = t29 * qJD(1);
t494 = t336 * t337;
t31 = t142 * t494 + t148 * t497 - t152 * t490 - t98;
t501 = t31 * qJD(1);
t32 = -t162 * t490 + (t526 + (t142 - t161) * t335) * t338;
t500 = t32 * qJD(1);
t33 = (t138 - t517) * t338 + (t148 + t515) * t336;
t499 = t33 * qJD(1);
t487 = t338 * t230;
t486 = t338 * t256;
t484 = t580 * t338;
t34 = -t142 * t338 - t152 * t336 + (t169 * t338 - t170 * t336) * t335;
t483 = t34 * qJD(1);
t53 = -t161 * t338 + (-t392 + 0.2e1 * t266) * t336;
t482 = t53 * qJD(1);
t54 = t162 * t338 + (-t267 + t242) * t336;
t481 = t54 * qJD(1);
t461 = qJD(5) * t335;
t301 = t338 * t461;
t160 = (-0.1e1 + t281) * t338 * t336;
t457 = t160 * qJD(2);
t480 = t457 + t301;
t203 = t337 * t441 + t564 + t571;
t453 = t203 * qJD(2);
t479 = -t162 * qJD(4) - t453;
t176 = t199 * qJD(2);
t321 = t336 * qJD(5);
t478 = t176 + t321;
t229 = t491 / 0.2e1;
t477 = t229 - t266 / 0.2e1;
t247 = t333 * t334 + t583;
t469 = qJD(2) * t336;
t291 = t335 * t469;
t476 = t247 * qJD(5) + t291;
t200 = t335 * t391 + t568;
t178 = t200 * qJD(2);
t473 = t321 - t178;
t472 = qJD(1) * qJ(2);
t280 = -t334 + t583;
t240 = t280 * t335;
t471 = qJD(1) * t240;
t241 = t337 * t583 - t313;
t470 = qJD(1) * t241;
t468 = qJD(3) * t335;
t467 = qJD(3) * t337;
t466 = qJD(4) * t161;
t465 = qJD(4) * t257;
t464 = qJD(4) * t335;
t460 = qJD(6) * t335;
t459 = qJD(6) * t337;
t458 = qJD(6) * t338;
t454 = t202 * qJD(2);
t237 = t281 * t334;
t452 = t237 * qJD(1);
t239 = t281 * t338;
t451 = t239 * qJD(1);
t450 = t239 * qJD(3);
t449 = t280 * qJD(1);
t448 = t281 * qJD(3);
t447 = t336 * qJD(1);
t446 = t336 * qJD(3);
t442 = t338 * qJD(4);
t440 = -pkin(4) / 0.2e1 - pkin(5) / 0.2e1;
t264 = t307 / 0.2e1;
t439 = t264 - t309 / 0.2e1;
t437 = pkin(8) * t464;
t436 = pkin(8) * t462;
t434 = pkin(5) / 0.2e1 + t560;
t405 = -t490 / 0.2e1;
t406 = t494 / 0.2e1;
t429 = -t516 / 0.2e1 + t245 * t405 + pkin(8) * t406;
t227 = t242 / 0.2e1;
t428 = t227 + t585;
t427 = t326 / 0.2e1 + t477;
t426 = qJ(2) * t447;
t425 = qJ(2) * t444;
t422 = t335 * t442;
t421 = t337 * t442;
t420 = t337 * t458;
t419 = t335 * t462;
t418 = t336 * t460;
t294 = t335 * t467;
t293 = t335 * t445;
t417 = t336 * t443;
t416 = t336 * t444;
t411 = t498 / 0.2e1;
t410 = t497 / 0.2e1;
t407 = -t494 / 0.2e1;
t404 = t490 / 0.2e1;
t403 = t487 / 0.2e1;
t402 = -t486 / 0.2e1;
t401 = t486 / 0.2e1;
t400 = t335 * t560;
t398 = t484 / 0.2e1;
t395 = t569 - t333 / 0.2e1;
t197 = (-0.1e1 / 0.2e1 + t395) * t338;
t132 = qJD(1) * t197 - t294;
t226 = t395 * t338;
t390 = qJD(1) * t226 - t294;
t248 = t335 * qJD(1) * t313;
t389 = qJD(3) * t226 + t248;
t292 = t335 * t447;
t388 = qJD(4) * t199 + t292;
t387 = qJD(4) * t200 - t292;
t308 = qJD(4) + t447;
t386 = t335 * t300;
t384 = t582 + t434;
t382 = t581 - t507 / 0.2e1;
t381 = t505 / 0.2e1 - t73 / 0.2e1;
t378 = 0.2e1 * t386;
t376 = -t487 / 0.2e1 + t579;
t375 = -qJ(5) * t463 - t454;
t373 = pkin(3) * t410 + pkin(8) * t411;
t371 = t509 + t558;
t21 = t138 * t148 + t142 * t152 + t169 * t170;
t370 = t21 * qJD(1) + t14 * qJD(2);
t348 = t578 - t356;
t16 = (-t510 / 0.2e1 - t434) * t338 + (t573 - t348) * t337 + t346;
t82 = -t220 * t335 + t230 * t337;
t369 = qJD(1) * t16 + qJD(3) * t82;
t23 = t508 / 0.2e1 + t376 * t337 + t348 * t335 + t428;
t81 = t220 * t337 + t510;
t368 = -qJD(1) * t23 + qJD(3) * t81;
t367 = t521 + t523;
t144 = t255 * t335 + t257 * t337;
t139 = t245 * t337 + t256 * t335;
t228 = -t242 / 0.2e1;
t344 = pkin(8) * t412 + (t245 * t561 + t576) * t335 + (t402 - t170 / 0.2e1) * t337;
t36 = t228 + t344 - t585;
t366 = -qJD(1) * t36 + qJD(3) * t139;
t140 = -t245 * t335 + t256 * t337;
t276 = pkin(8) * t407;
t353 = t335 * t401 + t516 / 0.2e1 + t245 * t404 + t276;
t38 = -t326 + t572 + (t576 + t573) * t337 + t353;
t365 = -qJD(1) * t38 + qJD(3) * t140;
t102 = (t337 * t384 + t323) * t338;
t165 = -t335 * t384 + t492;
t364 = qJD(1) * t102 - qJD(3) * t165;
t362 = -t148 * qJ(5) / 0.2e1 + t152 * t582;
t361 = t79 * t559 + t96 * qJ(5) / 0.2e1;
t120 = t227 + t373;
t360 = pkin(3) * t467 - qJD(1) * t120;
t121 = t276 + (-t555 / 0.2e1 + t573) * t337;
t359 = pkin(3) * t468 - qJD(1) * t121;
t358 = t220 * t578 + t230 * t579;
t46 = t427 + t429;
t355 = -qJD(1) * t46 + t245 * t468;
t192 = t308 * t490;
t234 = t335 * t444 - t467;
t238 = t282 * t334;
t352 = qJD(1) * t238 + t378;
t350 = -qJD(2) * t239 - t301 * t336;
t1 = t255 * t430 - t257 * t431 + t358 + t361;
t50 = t220 * t230;
t52 = -t307 / 0.2e1 + (t230 / 0.2e1 + t400) * t338;
t349 = -t1 * qJD(1) + t52 * qJD(2) + t50 * qJD(3);
t196 = (-0.1e1 / 0.2e1 + t333 / 0.2e1 + t569) * t336;
t25 = t493 / 0.2e1 + (t409 + t382) * t337 + (t336 * t440 + t381) * t335;
t347 = qJD(1) * t25 - qJD(2) * t196 - qJD(3) * t144;
t117 = t401 + t439;
t341 = (t379 - t588) * pkin(8) + t256 * t577 + t217 * t245 / 0.2e1;
t6 = t341 + t362;
t345 = t245 * t256 * qJD(3) + t6 * qJD(1) - t117 * qJD(2);
t126 = qJD(1) * t247 + t386 + t463;
t331 = qJ(5) * qJD(5);
t320 = t332 * qJD(5);
t312 = qJ(5) * t321;
t311 = t323 / 0.2e1;
t310 = t443 / 0.2e1;
t299 = t337 * t447;
t298 = t337 * t469;
t297 = t336 * t445;
t295 = t335 * t446;
t249 = t337 * t301;
t246 = t308 * qJ(5);
t236 = t424 + t468;
t235 = t299 + t462;
t233 = (t447 + qJD(4) / 0.2e1) * t338;
t216 = -t295 + t421;
t215 = -t337 * t446 - t422;
t214 = qJD(3) * t332 + t385;
t213 = t226 * qJD(4);
t204 = t335 * t440 + t400;
t201 = (0.1e1 + t281) * t566;
t198 = t332 * t562 + t333 * t561 + t562;
t191 = t337 * t416 + t295;
t190 = t308 * t497;
t189 = t234 * t336;
t183 = t203 * qJD(5);
t180 = t202 * qJD(5);
t175 = t198 * qJD(5);
t174 = t197 * qJD(5);
t158 = t337 * t398 + t405 * t580;
t149 = t160 * qJD(3);
t131 = -qJD(4) * t202 - t299;
t130 = -qJD(4) * t203 + t299;
t118 = t402 + t439;
t112 = qJD(3) * t198 - t248;
t111 = qJD(3) * t197 + t248;
t76 = pkin(3) * t405 + t229 - t266 + t276;
t75 = -t267 + t228 + t373;
t51 = -t410 * t580 + t264 + t403;
t45 = -t491 / 0.2e1 + t429 + t438;
t40 = t220 * t404 + t435 + (t573 + t408) * t337 - t357 + t438;
t37 = -t512 / 0.2e1 + t326 + t353 + t477;
t35 = t344 + t428;
t30 = t161 * t564 - t527 / 0.2e1 + t415 + t524 / 0.2e1 + t260 + pkin(4) * t406;
t26 = -t493 / 0.2e1 + qJ(5) * t407 + t382 * t337 + t381 * t335 + t580 * t411;
t24 = -t508 / 0.2e1 + t337 * t403 + t335 * t414 + t518 / 0.2e1 + t528 / 0.2e1 + qJ(6) * t412 + t428;
t19 = t311 + t556 / 0.2e1 + t343;
t17 = t413 + t398 + qJ(6) * t407 + t554 / 0.2e1 + (t414 + t166 / 0.2e1) * t337 + t376 * t335 + t427;
t12 = t14 * qJD(3);
t10 = t311 + t489 / 0.2e1 + t342;
t9 = t115 * t564 + t547 / 0.2e1 + t116 * t568 - t550 / 0.2e1 + t261 + t337 * t399;
t5 = t341 - t362;
t3 = t257 * t567 + t255 * t564 + (t92 * t563 + t73 * t567 + t136 / 0.2e1) * t338 + (t563 * t96 + t567 * t79 + t579) * t336;
t2 = t116 * t575 + t255 * t581 - t358 + t361 + t589 * t257 / 0.2e1;
t57 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t417, t280 * qJD(3), 0, 0, 0, qJ(2) * t443 + t469, -qJ(2) * t446 + qJD(2) * t338, -t333 * t417 - t334 * t419, -t238 * qJD(4) + t336 * t378, -qJD(3) * t241 - t336 * t422, qJD(3) * t240 - t336 * t421, t417, qJD(3) * t53 + qJD(4) * t84 + t298, -qJD(3) * t54 - qJD(4) * t83 - t291, qJD(3) * t34 + qJD(4) * t47 - t293 * t334 + t298, -qJD(3) * t31 - qJD(4) * t32 + t350, qJD(3) * t33 + qJD(4) * t48 + t476, qJD(2) * t55 + qJD(3) * t21 + qJD(4) * t28 + qJD(5) * t56, -qJD(3) * t22 - qJD(4) * t42 + t298 + (-t334 * t461 + t336 * t458) * t337, qJD(3) * t27 + qJD(4) * t41 + t338 * t418 + t476, qJD(3) * t15 + qJD(4) * t18 + qJD(6) * t237 - t350, qJD(2) * t44 + qJD(3) * t7 + qJD(4) * t13 + qJD(5) * t49 + qJD(6) * t43; 0, 0, 0, 0, qJD(1), t472, 0, 0, 0, 0, 0, t447, t444, 0, 0, 0, 0, 0, t130, t387, t130, -t451, -t387, qJD(4) * t19 + t12 + t183 + t537, t130, -t387, t451, qJD(3) * t3 + qJD(4) * t10 + t183 + t541; 0, 0, 0, 0, 0, 0, -t416, t449, -t446, -t443, 0, -t340 * t446 + t425, -t340 * t443 - t426, -t213 + (-t333 * t444 - t294) * t336, t336 * t351 - 0.2e1 * t338 * t419, t296 - t470, t300 + t471, t233, t482 + (t335 * t380 - t433) * qJD(3) + t76 * qJD(4), -t481 + (-pkin(8) * t490 + (pkin(3) * t337 + t496) * t336) * qJD(3) + t75 * qJD(4), t483 + (-t335 * t371 - t517) * qJD(3) + t37 * qJD(4) + t175, qJD(3) * t367 + t30 * qJD(4) - t501, t499 + (-t169 * t335 + t337 * t371) * qJD(3) + t35 * qJD(4) + t249 (pkin(8) * t367 + t169 * t245) * qJD(3) + t5 * qJD(4) + t45 * qJD(5) + t370, -t511 + (t220 * t498 - t507 + t530) * qJD(3) + t17 * qJD(4) + t175 + t418, t504 + (-t220 * t494 + t505 + t531) * qJD(3) + t24 * qJD(4) + t249 - t336 * t459, t522 + ((-t96 + t508) * t337 + (-t79 - t506) * t335) * qJD(3) + t9 * qJD(4), t551 + t3 * qJD(2) + (t136 * t220 + t255 * t79 + t257 * t96) * qJD(3) + t2 * qJD(4) + t40 * qJD(5) + t26 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t389, -t352, -t190, -t192, t310, qJD(3) * t76 + t479 + t534, qJD(3) * t75 + t178 + t466 - t535, qJD(3) * t37 + t479 + t540, t30 * qJD(3) - qJD(4) * t474 - t301 - t500, qJD(3) * t35 - t466 + t473 + t539, t503 + t19 * qJD(2) + t5 * qJD(3) + (-pkin(4) * t162 - qJ(5) * t161) * qJD(4) + t138 * qJD(5), qJD(3) * t17 - qJD(4) * t116 - t453 - t543, qJD(3) * t24 + qJD(4) * t115 + t473 + t544, t514 + t9 * qJD(3) + (-t335 * t484 + t307) * qJD(4) + t301, t532 + t10 * qJD(2) + t2 * qJD(3) + (qJ(5) * t115 - t116 * t580) * qJD(4) + t92 * qJD(5) + t158 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t190, t126, qJD(3) * t45 + qJD(4) * t138 + t453 + t536, t112, t126, t190, qJD(3) * t40 + qJD(4) * t92 + t453 + t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, t189, t452, qJD(3) * t26 + qJD(4) * t158 + t542; 0, 0, 0, 0, -qJD(1), -t472, 0, 0, 0, 0, 0, -t447, -t444, 0, 0, 0, 0, 0, t131, t388, t131, t451, -t388, qJD(4) * t20 + t12 + t180 - t537, t131, -t388, -t451, qJD(3) * t4 + qJD(4) * t11 + t180 - t541; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, 0, 0, 0, t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t446, -t443, 0, 0, 0, 0, 0, t215, -t216, t215, t450, t216, t525 + (pkin(8) * t239 + t509) * qJD(3) + t118 * qJD(4) + t480, t215, t216, -t450, t552 + (t144 * t338 - t220 * t336) * qJD(3) + t51 * qJD(4) + t201 * qJD(6) + t480; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t587, -t86, 0, -t587, t118 * qJD(3) - t336 * t354 + t297 + t513, -t86, -t587, 0, t51 * qJD(3) + t463 * t590 + t297 + t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201 * qJD(3); 0, 0, 0, 0, 0, 0, t416, -t449, 0, 0, 0, -t425, t426, t333 * t416 - t213, -0.2e1 * t335 * t192, t302 + t470, -t423 - t471, -t233, qJD(4) * t121 - t482, qJD(4) * t120 + t481, qJD(4) * t38 - t174 - t483, -qJD(4) * t29 + t297 + t501, qJD(4) * t36 + t249 - t499, qJD(4) * t6 + qJD(5) * t46 - t370, qJD(4) * t16 - t174 + t511, -qJD(4) * t23 + t249 - t504, -qJD(4) * t8 - t297 - t522, -qJD(2) * t4 - qJD(4) * t1 - qJD(5) * t39 + qJD(6) * t25 - t551; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t117 - t457 - t525, 0, 0, 0, qJD(4) * t52 - qJD(6) * t196 - t457 - t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419, t282 * qJD(4), 0, 0, 0, -pkin(3) * t464, -pkin(3) * t462, -qJD(4) * t140 + t293, 0, -qJD(4) * t139 + t320 (qJD(4) * t256 - t461) * t245, qJD(4) * t82 + t293, qJD(4) * t81 + t320, qJD(6) * t281, qJD(4) * t50 - qJD(6) * t144 + t220 * t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t390, -t351, t235, -t292 - t464, -t444 / 0.2e1, -t359 - t436, -t360 + t437, -t365 - t436, -t502 + t591, -t366 - t437, pkin(8) * t591 + t345, t369 - t465, -qJD(4) * t255 + t368, -qJD(4) * t590 - t445 - t548 (-qJ(5) * t255 - t257 * t580) * qJD(4) + t257 * qJD(5) + t204 * qJD(6) + t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, t235, t214, -t355 + t436, -t132, t214, -t235, t220 * t468 + t465 - t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t448, qJD(4) * t204 + t347; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t389, t352, t189, t191, t310, -qJD(3) * t121 + t454 - t534, -qJD(3) * t120 - t176 + t535, -qJD(3) * t38 + t454 - t540, qJD(3) * t29 + t500, -qJD(3) * t36 + t478 - t539, -qJD(2) * t20 - qJD(3) * t6 + t312 - t503, -qJD(3) * t16 + t420 + t454 + t543, qJD(3) * t23 + t335 * t458 + t478 - t544, qJD(3) * t8 - t514, -qJD(2) * t11 + qJD(3) * t1 + qJD(6) * t102 + t312 - t532; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t455, -t456, t455, 0, t456, qJD(3) * t117 - t513, t455, t456, 0, -qJD(3) * t52 - t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t390, t351, -t299, t292, t444 / 0.2e1, t359, t360, t365, t502, t366, -t345, -t369 + t460, -t368 - t459, t548, -qJD(6) * t165 - t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t331, 0, qJD(5), 0, t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t308, t246, 0, t308, 0, t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, t234, 0, t364; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t189, -t126, -qJD(3) * t46 + t375 - t536, t111, -t126, -t189, qJD(3) * t39 + t375 - t420 - t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t455, 0, 0, 0, -t455; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, -t299, -t214, t355, t132, -t214, t299, t545 + (-qJD(3) * t220 - qJD(6)) * t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t308, -t246, 0, -t308, 0, -t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, -t190, -t452, -qJD(3) * t25 - qJD(4) * t102 + t338 * t445 - t542; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t464, t462, -t448, qJD(4) * t165 - t347 + t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, -t234, 0, -t364; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t57;
