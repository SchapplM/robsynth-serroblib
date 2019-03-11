% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x29]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:04
% EndTime: 2019-03-09 05:21:18
% DurationCPUTime: 7.20s
% Computational Cost: add. (9834->384), mult. (17491->522), div. (0->0), fcn. (20828->8), ass. (0->324)
t301 = sin(qJ(4));
t302 = sin(qJ(3));
t304 = cos(qJ(4));
t305 = cos(qJ(3));
t276 = t301 * t302 - t304 * t305;
t277 = t301 * t305 + t302 * t304;
t299 = sin(pkin(10));
t477 = cos(pkin(10));
t246 = -t299 * t276 + t477 * t277;
t535 = t477 * t246;
t266 = t477 * t276;
t445 = t299 * t277;
t333 = t266 + t445;
t545 = t333 * t299;
t570 = (t545 / 0.2e1 + t535 / 0.2e1) * pkin(4);
t487 = pkin(4) * qJD(4);
t569 = t487 * (t535 + t545);
t292 = pkin(3) * t304 + pkin(4);
t444 = t299 * t301;
t264 = -pkin(3) * t444 + t292 * t477;
t536 = t264 * t246;
t359 = t477 * t301;
t265 = pkin(3) * t359 + t299 * t292;
t546 = t333 * t265;
t567 = qJD(3) * (t546 + t536);
t300 = sin(qJ(6));
t306 = -pkin(1) - pkin(7);
t531 = -pkin(8) + t306;
t281 = t531 * t302;
t282 = t531 * t305;
t201 = -t304 * t281 - t301 * t282;
t329 = qJ(5) * t277 + t201;
t200 = t301 * t281 - t304 * t282;
t528 = t276 * qJ(5) - t200;
t543 = t299 * t528 - t329 * t477;
t465 = t543 * t333;
t303 = cos(qJ(6));
t290 = pkin(3) * t302 + qJ(2);
t256 = pkin(4) * t277 + t290;
t315 = pkin(5) * t246 + pkin(9) * t333 + t256;
t564 = t543 * t300;
t55 = -t303 * t315 + t564;
t566 = -t300 * t465 + t55 * t333;
t563 = t543 * t303;
t56 = t300 * t315 + t563;
t565 = -t303 * t465 + t56 * t333;
t559 = t536 / 0.2e1 + t546 / 0.2e1;
t505 = -t333 / 0.2e1;
t557 = t505 * t543;
t492 = t303 / 0.2e1;
t493 = -t303 / 0.2e1;
t542 = t299 * t329 + t528 * t477;
t556 = (-t492 + t493) * t542;
t494 = t300 / 0.2e1;
t495 = -t300 / 0.2e1;
t555 = (-t494 + t495) * t542;
t552 = t333 ^ 2;
t551 = t246 / 0.2e1;
t549 = t246 * t333;
t274 = (t304 * t477 - t444) * pkin(3);
t499 = -t274 / 0.2e1;
t548 = t246 * t499;
t547 = t300 * t333;
t440 = t303 * t333;
t538 = t333 / 0.2e1;
t544 = t333 * t538;
t534 = qJD(1) * t246;
t541 = -qJD(6) - t534;
t388 = qJD(3) + qJD(4);
t540 = -0.2e1 * t333;
t507 = -t246 / 0.2e1;
t357 = t246 * t388;
t533 = -0.2e1 * t300 * t440;
t527 = t246 * t551 + t544;
t491 = t276 * pkin(4);
t120 = -pkin(5) * t333 + pkin(9) * t246 - t491;
t295 = t305 * pkin(3);
t119 = t120 + t295;
t261 = pkin(9) + t265;
t260 = -pkin(5) - t264;
t504 = -t260 / 0.2e1;
t339 = t261 * t551 - t333 * t504;
t325 = t119 / 0.2e1 + t339;
t23 = t303 * t325;
t242 = -t445 / 0.2e1 - t266 / 0.2e1;
t316 = t538 + t242;
t141 = t316 * t300;
t407 = t141 * qJD(2);
t417 = qJD(3) * t300;
t526 = qJD(1) * t23 - t260 * t417 + t407;
t21 = t300 * t325;
t146 = t316 * t303;
t405 = t146 * qJD(2);
t416 = qJD(3) * t303;
t525 = -qJD(1) * t21 - t260 * t416 + t405;
t273 = (t299 * t304 + t359) * pkin(3);
t338 = t273 * t505 + t548;
t453 = t246 * t260;
t288 = t299 * pkin(4) + pkin(9);
t458 = t333 * t288;
t459 = t333 * t261;
t524 = t458 / 0.2e1 + t459 / 0.2e1 - t453 / 0.2e1 + t338;
t523 = t388 * t201;
t522 = t388 * t200;
t289 = -pkin(4) * t477 - pkin(5);
t497 = -t289 / 0.2e1;
t351 = t499 + t497 + t504;
t190 = t351 * t303;
t498 = t288 / 0.2e1;
t337 = t246 * t498 - t333 * t497;
t324 = t120 / 0.2e1 + t337;
t25 = t300 * t324;
t413 = qJD(4) * t303;
t521 = -qJD(1) * t25 + qJD(3) * t190 - t289 * t413 + t405;
t189 = t351 * t300;
t27 = t303 * t324;
t414 = qJD(4) * t300;
t520 = qJD(1) * t27 + qJD(3) * t189 - t289 * t414 + t407;
t11 = t538 * t543 + t557;
t519 = t11 * qJD(1);
t512 = t543 / 0.2e1;
t10 = t333 * t512 + t557;
t518 = t10 * qJD(1);
t298 = t303 ^ 2;
t514 = t300 ^ 2;
t136 = (t514 / 0.2e1 - t298 / 0.2e1) * t333;
t419 = qJD(1) * t303;
t378 = t300 * t419;
t515 = -t136 * t388 + t378 * t552;
t386 = -t298 + t514;
t169 = t378 * t540 + t386 * t388;
t513 = t11 * qJD(3) + t10 * qJD(4);
t503 = t264 / 0.2e1;
t502 = -t265 / 0.2e1;
t500 = -t273 / 0.2e1;
t496 = t299 / 0.2e1;
t313 = -t549 / 0.2e1 + t333 * t551;
t39 = t313 * t303;
t490 = t388 * t39;
t66 = t313 * t300;
t489 = t388 * t66;
t488 = pkin(3) * qJD(4);
t486 = qJD(3) * pkin(3);
t441 = t303 * t246;
t1 = t119 * t441 + t566;
t485 = t1 * qJD(1);
t443 = t300 * t246;
t2 = -t119 * t443 + t565;
t484 = t2 * qJD(1);
t3 = t120 * t441 + t566;
t483 = t3 * qJD(1);
t4 = -t120 * t443 + t565;
t482 = t4 * qJD(1);
t438 = t39 * qJD(2);
t478 = -qJD(5) * t547 + t438;
t30 = t246 * t55 - t542 * t547;
t476 = qJD(1) * t30;
t31 = -t246 * t56 + t440 * t542;
t475 = qJD(1) * t31;
t50 = -t246 * t543 + t333 * t542;
t474 = qJD(1) * t50;
t456 = t246 ^ 2;
t346 = -t456 + t552;
t67 = t346 * t300;
t473 = qJD(1) * t67;
t76 = t456 + t552;
t68 = t76 * t300;
t472 = qJD(1) * t68;
t308 = -(t260 / 0.2e1 + t497) * t246 - (-t261 / 0.2e1 + t498) * t333 + t338;
t366 = t512 - t543 / 0.2e1;
t14 = t300 * t308 - t303 * t366;
t461 = t14 * qJD(1);
t448 = t298 * t333;
t439 = t39 * qJD(1);
t311 = t338 + t559;
t40 = t570 - t311;
t437 = t40 * qJD(1);
t48 = t549 / 0.2e1 + t333 * t507;
t436 = t48 * qJD(1);
t435 = t48 * qJD(2);
t309 = 0.1e1 / 0.2e1 + t527;
t58 = t309 * t300;
t434 = t58 * qJD(1);
t59 = t309 * t303;
t433 = t59 * qJD(1);
t432 = t66 * qJD(1);
t431 = t66 * qJD(2);
t69 = t346 * t303;
t430 = t69 * qJD(1);
t70 = t76 * t303;
t429 = t70 * qJD(1);
t314 = t246 * t507 - t544;
t72 = -0.1e1 / 0.2e1 + t314;
t428 = t72 * qJD(1);
t426 = t76 * qJD(1);
t340 = -t246 * t502 - t333 * t503;
t268 = -t491 / 0.2e1;
t387 = t268 + t295 / 0.2e1;
t83 = t340 + t387;
t425 = t83 * qJD(1);
t424 = (t413 + t416) * t333;
t423 = qJD(1) * qJ(2);
t420 = qJD(1) * t256;
t418 = qJD(2) * t246;
t415 = qJD(4) * t290;
t412 = qJD(6) * t300;
t411 = qJD(6) * t303;
t384 = -t477 / 0.2e1;
t331 = -t246 * t496 - t333 * t384;
t115 = (t276 / 0.2e1 + t331) * pkin(4);
t410 = t115 * qJD(1);
t409 = t136 * qJD(1);
t408 = t547 * qJD(1);
t363 = 0.2e1 * t551;
t139 = t363 * t300;
t130 = t139 * qJD(1);
t143 = t363 * t303;
t132 = t143 * qJD(1);
t231 = t246 * t492;
t367 = -t441 / 0.2e1;
t144 = t367 - t231;
t406 = t144 * qJD(1);
t151 = t386 * t552;
t404 = t151 * qJD(1);
t215 = -t276 ^ 2 + t277 ^ 2;
t399 = t215 * qJD(1);
t238 = -t276 * t290 + t277 * t295;
t398 = t238 * qJD(1);
t239 = -t276 * t295 - t277 * t290;
t397 = t239 * qJD(1);
t396 = t242 * qJD(1);
t395 = t276 * qJD(1);
t394 = t277 * qJD(1);
t285 = t302 ^ 2 - t305 ^ 2;
t393 = t285 * qJD(1);
t392 = t302 * qJD(1);
t391 = t302 * qJD(3);
t390 = t305 * qJD(1);
t389 = t305 * qJD(3);
t383 = qJ(2) * t392;
t382 = qJ(2) * t390;
t381 = t300 * t534;
t380 = t246 * t419;
t379 = qJD(1) * t448;
t377 = qJD(6) * t549;
t376 = t290 * t395;
t375 = t290 * t394;
t287 = t300 * t411;
t374 = t302 * t390;
t369 = t289 * t494;
t358 = pkin(3) * t388;
t356 = t388 * t303;
t251 = t388 * t276;
t355 = t388 * t300;
t352 = qJD(5) * t440 + t431;
t350 = t300 * t356;
t349 = t303 * t355;
t20 = t256 * t491;
t348 = -t20 * qJD(1) + t10 * qJD(2);
t19 = t256 * (t295 - t491);
t347 = t19 * qJD(1) + t11 * qJD(2);
t345 = -t453 + t459;
t344 = -t246 * t289 + t458;
t343 = t333 * t541;
t17 = t300 * t366 + t303 * t308;
t336 = -qJD(1) * t17 - t273 * t417;
t90 = qJD(6) * t242 - t333 * t534;
t332 = t350 * t540;
t168 = -t264 * t273 + t265 * t274;
t307 = t333 * t500 + t548 + t559;
t41 = -t570 + t307;
t312 = (t499 + t503) * t543 + (-t500 + t502) * t542;
t323 = (t384 * t543 + t496 * t542) * pkin(4);
t5 = t323 + t312;
t328 = t5 * qJD(1) + t41 * qJD(2) - t168 * qJD(3);
t317 = t333 * t492;
t284 = t386 * qJD(6);
t259 = t273 * t414;
t253 = t276 * t394;
t252 = t388 * t277;
t237 = t333 * t419;
t192 = t274 * t493 + (t260 + t289) * t492;
t191 = t260 * t494 + t274 * t495 + t369;
t149 = t388 * t242;
t148 = t317 + t440 / 0.2e1;
t147 = t333 * t494 + t547 / 0.2e1;
t145 = (t507 + t551) * t303;
t142 = t441 / 0.2e1 - t231;
t140 = -t443 / 0.2e1 - t246 * t495;
t129 = t140 * qJD(6);
t128 = t139 * qJD(6);
t123 = t136 * qJD(6);
t122 = t132 + t411;
t121 = -t130 - t412;
t116 = pkin(4) * t331 + t268;
t114 = t349 + t409;
t113 = -t350 - t409;
t85 = t541 * t533;
t84 = -t340 + t387;
t71 = 0.1e1 / 0.2e1 + t314;
t61 = t246 * t367 - t317 * t333 + t492;
t60 = t300 * t527 + t495;
t57 = -t246 * t379 + t123;
t54 = qJD(6) * t143 - t430;
t53 = -t128 + t473;
t52 = -qJD(6) * t141 + t432;
t51 = t123 - (t349 - t379) * t246;
t49 = (-qJD(6) + t534) * t533 + t386 * t357;
t43 = -t570 - t307;
t42 = t570 + t311;
t34 = t142 * qJD(6) - t333 * t355 + t430;
t33 = t129 - t424 - t473;
t32 = t147 * qJD(6) - t246 * t356 - t432;
t29 = -qJD(6) * t146 + t439;
t28 = t120 * t492 - t303 * t337 + t555;
t26 = t120 * t495 + t300 * t337 + t556;
t24 = t119 * t492 - t303 * t339 + t555;
t22 = t119 * t495 + t300 * t339 + t556;
t18 = t148 * qJD(6) + t246 * t355 - t439;
t16 = -t289 * t231 + t524 * t303 + t564;
t15 = -t246 * t369 + t524 * t300 - t563;
t7 = t388 * t48;
t6 = t323 - t312;
t8 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t302 * t389, t285 * qJD(3), 0, 0, 0, qJ(2) * t389 + qJD(2) * t302, -qJ(2) * t391 + qJD(2) * t305, t277 * t251, t388 * t215, 0, 0, 0, qJD(2) * t277 + qJD(3) * t238 - t276 * t415, -qJD(2) * t276 + qJD(3) * t239 - t277 * t415, qJD(5) * t76, qJD(2) * t256 + qJD(3) * t19 - qJD(4) * t20 + qJD(5) * t50, -t287 * t552 + t357 * t448, t151 * qJD(6) + t246 * t332, t300 * t377 + t388 * t69, t303 * t377 - t388 * t67, -t333 * t357, qJD(3) * t1 + qJD(4) * t3 + qJD(5) * t68 + qJD(6) * t31 + t303 * t418, qJD(3) * t2 + qJD(4) * t4 + qJD(5) * t70 + qJD(6) * t30 - t300 * t418; 0, 0, 0, 0, qJD(1), t423, 0, 0, 0, 0, 0, t392, t390, 0, 0, 0, 0, 0, t394, -t395, t7, qJD(5) * t71 + t420 + t513, 0, 0, 0, 0, 0, qJD(6) * t61 + t380 - t489, qJD(6) * t60 - t381 - t490; 0, 0, 0, 0, 0, 0, -t374, t393, -t391, -t389, 0, -t306 * t391 + t382, -t306 * t389 - t383, t253, t399, -t252, t251, 0, t398 + t523, t397 + t522, t42 * qJD(4) + t435 + t567 (-t264 * t543 + t265 * t542) * qJD(3) + t6 * qJD(4) + t84 * qJD(5) + t347, t51, t49, t34, t33, t90, t485 - t431 + (t300 * t345 - t563) * qJD(3) + t15 * qJD(4) + t24 * qJD(6), t484 + (t303 * t345 + t564) * qJD(3) + t16 * qJD(4) + t22 * qJD(6) - t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t399, -t252, t251, 0, -t376 + t523, -t375 + t522, t42 * qJD(3) + t435 + t569, t6 * qJD(3) + (t299 * t542 - t477 * t543) * t487 + t116 * qJD(5) + t348, t51, t49, t34, t33, t90, t483 - t431 + t15 * qJD(3) + (t300 * t344 - t563) * qJD(4) + t28 * qJD(6), t482 + t16 * qJD(3) + (t303 * t344 + t564) * qJD(4) + t26 * qJD(6) - t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t426, qJD(2) * t71 + qJD(3) * t84 + qJD(4) * t116 + t474, 0, 0, 0, 0, 0, t129 + t472, t145 * qJD(6) + t429; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t515, -t332 + t404, t142 * t388 - t300 * t343, t140 * t388 - t303 * t343, t149, qJD(2) * t61 + qJD(3) * t24 + qJD(4) * t28 + qJD(5) * t140 - qJD(6) * t56 + t475, qJD(2) * t60 + qJD(3) * t22 + qJD(4) * t26 + qJD(5) * t145 + qJD(6) * t55 + t476; 0, 0, 0, 0, -qJD(1), -t423, 0, 0, 0, 0, 0, -t392, -t390, 0, 0, 0, 0, 0, -t394, t395, t7, qJD(5) * t72 - t420 + t513, 0, 0, 0, 0, 0, -qJD(6) * t59 - t380 - t489, qJD(6) * t58 + t381 - t490; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391, -t389, 0, 0, 0, 0, 0, -t252, t251, t436, t43 * qJD(4) + t519 - t567, 0, 0, 0, 0, 0, t32, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t252, t251, t436, t43 * qJD(3) + t518 - t569, 0, 0, 0, 0, 0, t32, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t428, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147 * t388 - t246 * t411 - t433, t148 * t388 + t246 * t412 + t434; 0, 0, 0, 0, 0, 0, t374, -t393, 0, 0, 0, -t382, t383, -t253, -t399, 0, 0, 0, -t398, -t397, -qJD(4) * t40 - t435, -qJD(4) * t5 - qJD(5) * t83 - t347, t57, t85, t54, t53, -t90, qJD(4) * t14 - qJD(6) * t23 + t352 - t485, qJD(4) * t17 + qJD(6) * t21 + t478 - t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t436, -qJD(4) * t41 - t519, 0, 0, 0, 0, 0, t52, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301 * t488, -t304 * t488, 0, t168 * qJD(4), t287, -t284, 0, 0, 0, t260 * t412 - t273 * t413, t260 * t411 + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301 * t358, -t304 * t358, -t437 (-t273 * t477 + t274 * t299) * t487 - t328, t287, -t284, 0, 0, 0, t191 * qJD(6) - t273 * t356 + t461, qJD(6) * t192 + t259 - t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t425, 0, 0, 0, 0, 0, t237, -t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t169, t122, t121, -t396, qJD(4) * t191 - t261 * t411 - t526, qJD(4) * t192 + t261 * t412 - t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253, -t399, 0, 0, 0, t376, t375, qJD(3) * t40 - t435, qJD(3) * t5 + qJD(5) * t115 - t348, t57, t85, t54, t53, -t90, -qJD(3) * t14 - qJD(6) * t27 + t352 - t483, -qJD(3) * t17 + qJD(6) * t25 + t478 - t482; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t436, qJD(3) * t41 - t518, 0, 0, 0, 0, 0, t52, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301 * t486, t304 * t486, t437, t328, t287, -t284, 0, 0, 0, -qJD(6) * t189 + t273 * t416 - t461, -qJD(6) * t190 + t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, -t284, 0, 0, 0, t289 * t412, t289 * t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t410, 0, 0, 0, 0, 0, t237, -t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t169, t122, t121, -t396, -t288 * t411 - t520, t288 * t412 - t521; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, -qJD(2) * t72 + qJD(3) * t83 - qJD(4) * t115 - t474, 0, 0, 0, 0, 0, -t128 - t424 - t472, t144 * qJD(6) + t388 * t547 - t429; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t428, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t425, 0, 0, 0, 0, 0, -t237, t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t410, 0, 0, 0, 0, 0, -t237, t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t406 - t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t515, t349 * t540 - t404, -t143 * t388 - t333 * t381, t139 * t388 - t333 * t380, t149, qJD(2) * t59 + qJD(3) * t23 + qJD(4) * t27 + qJD(5) * t139 - t475, -qJD(2) * t58 - qJD(3) * t21 - qJD(4) * t25 - qJD(5) * t144 - t476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141 * t388 + t433, t146 * t388 - t434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t169, -t132, t130, t396, qJD(4) * t189 + t526, qJD(4) * t190 + t525; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t169, -t132, t130, t396, t520, t521; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t8;
