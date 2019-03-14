% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRPR5
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
% cmat_reg [(6*%NQJ)%x32]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:36
% EndTime: 2019-03-09 05:13:54
% DurationCPUTime: 9.04s
% Computational Cost: add. (8626->446), mult. (16982->545), div. (0->0), fcn. (20338->8), ass. (0->338)
t294 = sin(qJ(6));
t420 = qJD(6) * t294;
t292 = sin(pkin(10));
t293 = cos(pkin(10));
t296 = sin(qJ(3));
t298 = cos(qJ(3));
t254 = t292 * t298 + t293 * t296;
t506 = cos(qJ(4));
t250 = t506 * t254;
t253 = -t296 * t292 + t293 * t298;
t295 = sin(qJ(4));
t450 = t295 * t253;
t530 = t250 + t450;
t550 = t294 * t530;
t566 = -t550 / 0.2e1;
t575 = 0.2e1 * t566;
t577 = t575 * qJD(1);
t588 = t577 - t420;
t288 = qJD(3) + qJD(4);
t500 = pkin(7) + qJ(2);
t258 = t500 * t292;
t259 = t500 * t293;
t333 = t298 * t258 + t296 * t259;
t164 = -t254 * pkin(8) - t333;
t332 = t258 * t296 - t259 * t298;
t165 = pkin(8) * t253 - t332;
t95 = t506 * t164 - t295 * t165;
t70 = -pkin(5) * t530 + t95;
t587 = t70 * t294;
t297 = cos(qJ(6));
t586 = t70 * t297;
t558 = t530 ^ 2;
t249 = t506 * t253;
t449 = t295 * t254;
t359 = -t249 + t449;
t559 = t359 ^ 2;
t562 = -t559 + t558;
t572 = t562 * t297;
t585 = qJD(1) * t572;
t573 = t562 * t294;
t584 = qJD(1) * t573;
t548 = t297 * t359;
t568 = -t548 / 0.2e1;
t570 = t548 / 0.2e1 + t568;
t583 = qJD(3) * t570;
t582 = qJD(4) * t570;
t581 = qJD(5) * t570;
t565 = t550 / 0.2e1;
t571 = t566 + t565;
t580 = qJD(6) * t571;
t579 = qJD(6) * t575;
t409 = t570 * qJD(1);
t578 = t570 * qJD(2);
t574 = 0.2e1 * t568;
t576 = qJD(4) * t574 + t580;
t308 = t295 * t164;
t384 = t506 * t165;
t528 = t384 + t308;
t537 = -pkin(5) * t359 + t528;
t563 = t297 * t537;
t521 = t563 / 0.2e1;
t569 = t562 * qJD(1);
t564 = t95 * qJ(5);
t290 = t294 ^ 2;
t291 = t297 ^ 2;
t360 = t291 / 0.2e1 - t290 / 0.2e1;
t120 = t360 * t359;
t531 = t288 * t297;
t549 = t294 * t531;
t91 = -t120 * qJD(1) + t549;
t270 = t290 - t291;
t132 = t270 * t559;
t561 = t132 * qJD(1) + 0.2e1 * t359 * t549;
t560 = t288 * t95;
t501 = t295 * pkin(3);
t273 = qJ(5) + t501;
t514 = t273 / 0.2e1;
t518 = -t359 / 0.2e1;
t557 = t359 / 0.2e1;
t519 = t530 / 0.2e1;
t555 = t528 * pkin(4);
t554 = t530 * pkin(4);
t274 = -t293 * pkin(2) - pkin(1);
t242 = -t253 * pkin(3) + t274;
t552 = qJ(5) * t530;
t311 = t242 - t552;
t505 = pkin(4) * t359;
t92 = t311 + t505;
t553 = t359 * t92;
t496 = t530 * t92;
t388 = t506 * pkin(3);
t277 = -t388 - pkin(4);
t271 = -pkin(9) + t277;
t551 = t271 * t359;
t547 = t359 * qJ(5);
t546 = t359 * t514;
t520 = pkin(4) + pkin(9);
t446 = t359 * t520;
t544 = t520 * t530;
t462 = t530 * t273;
t305 = t308 / 0.2e1;
t43 = t305 - t308 / 0.2e1;
t443 = t43 * qJD(1);
t527 = qJD(3) * t501 + t443;
t350 = (t557 + t518) * t294;
t525 = t350 * qJD(1);
t543 = t531 + t525;
t542 = qJD(2) * t359;
t226 = t449 / 0.2e1 - t249 / 0.2e1;
t541 = qJD(6) * t226;
t540 = t226 * qJD(1);
t539 = t359 * qJD(1);
t538 = t530 * qJD(1);
t419 = qJD(6) * t297;
t373 = t294 * t419;
t173 = -0.2e1 * t359 * t373;
t429 = qJD(1) * t297;
t379 = t294 * t429;
t352 = t359 * t379;
t348 = -0.2e1 * t352;
t536 = t348 * t530 + t173;
t535 = t288 * t359;
t532 = t270 * t288;
t371 = -t552 / 0.2e1;
t529 = t446 / 0.2e1 + t371;
t499 = pkin(3) * qJD(4);
t280 = t295 * t499;
t526 = t280 + t527;
t524 = t350 * qJD(3);
t398 = t359 * qJD(5);
t509 = t295 / 0.2e1;
t304 = (t506 * t518 + t509 * t530) * pkin(3) - t462 / 0.2e1;
t302 = -t551 / 0.2e1 + t304;
t351 = t250 / 0.2e1;
t228 = t351 - t250 / 0.2e1;
t42 = t384 + 0.2e1 * t305;
t389 = -t228 * qJD(2) + t42 * qJD(3) + qJD(4) * t528;
t523 = t120 * t288 + t379 * t559;
t515 = -t530 / 0.2e1;
t513 = t277 / 0.2e1;
t510 = -t294 / 0.2e1;
t508 = -t297 / 0.2e1;
t507 = t297 / 0.2e1;
t502 = t254 * pkin(3);
t67 = t311 + t446;
t23 = t294 * t67 + t586;
t315 = t359 * t70;
t341 = t547 + t502;
t74 = t341 + t544;
t1 = t23 * t359 - t297 * t315 - t550 * t74;
t498 = t1 * qJD(1);
t24 = t297 * t67 - t587;
t447 = t297 * t530;
t2 = t24 * t359 + t294 * t315 - t447 * t74;
t497 = t2 * qJD(1);
t79 = t547 + t544;
t495 = t294 * t79;
t3 = -t530 * t495 + (t23 - t586) * t359;
t493 = t3 * qJD(1);
t4 = (t24 + t587) * t359 - t79 * t447;
t492 = t4 * qJD(1);
t100 = t341 + t554;
t7 = t100 * t92;
t489 = t7 * qJD(1);
t486 = t537 * t294;
t133 = t547 + t554;
t8 = t92 * t133;
t485 = t8 * qJD(1);
t146 = 0.2e1 * t351 + t450;
t484 = t146 * qJD(2) - t43 * qJD(3);
t19 = t23 * t530 + t537 * t548;
t480 = qJD(1) * t19;
t454 = t294 * t359;
t20 = -t24 * t530 + t454 * t537;
t479 = qJD(1) * t20;
t25 = -t359 * t528 - t530 * t95;
t478 = qJD(1) * t25;
t26 = -t100 * t359 - t496;
t477 = qJD(1) * t26;
t27 = -t100 * t530 + t553;
t476 = qJD(1) * t27;
t28 = -t133 * t359 - t496;
t475 = qJD(1) * t28;
t29 = -t133 * t530 + t553;
t474 = qJD(1) * t29;
t56 = t559 + t558;
t51 = t56 * t294;
t471 = qJD(1) * t51;
t53 = t56 * t297;
t469 = qJD(1) * t53;
t80 = -t242 * t530 - t359 * t502;
t466 = qJD(1) * t80;
t81 = t242 * t359 - t502 * t530;
t465 = qJD(1) * t81;
t452 = t294 * t297;
t300 = t277 * t518 + t304;
t323 = pkin(4) * t557 + t371;
t30 = t300 - t323;
t445 = t30 * qJD(1);
t361 = t514 + qJ(5) / 0.2e1;
t44 = -t502 / 0.2e1 - t361 * t359 + (t513 - pkin(4) / 0.2e1) * t530;
t442 = t44 * qJD(1);
t362 = 0.2e1 * t557;
t363 = 0.2e1 * t519;
t48 = pkin(4) * t363 + qJ(5) * t362;
t441 = t48 * qJD(1);
t439 = t56 * qJD(1);
t281 = qJD(4) * t388;
t284 = qJD(5) * t294;
t437 = t294 * t281 + t284;
t285 = qJD(5) * t297;
t436 = t297 * t281 + t285;
t282 = t294 * qJD(4);
t283 = t294 * qJD(3);
t434 = t282 + t283;
t268 = t292 ^ 2 + t293 ^ 2;
t102 = (t519 + t515) * t452;
t433 = qJD(1) * t102;
t431 = qJD(1) * t242;
t430 = qJD(1) * t294;
t427 = qJD(3) * t530;
t425 = qJD(3) * t297;
t423 = qJD(4) * t530;
t422 = qJD(4) * t242;
t421 = qJD(4) * t297;
t418 = qJD(6) * t520;
t110 = 0.2e1 * t565;
t417 = t110 * qJD(1);
t115 = t362 * t294;
t413 = t115 * qJD(1);
t121 = t363 * t297;
t107 = t121 * qJD(1);
t123 = t362 * t297;
t410 = t123 * qJD(1);
t408 = t548 * qJD(1);
t405 = t146 * qJD(1);
t172 = t253 ^ 2 - t254 ^ 2;
t404 = t172 * qJD(1);
t401 = t228 * qJD(1);
t211 = t228 * qJD(4);
t400 = t558 * qJD(1);
t396 = t253 * qJD(1);
t252 = t253 * qJD(3);
t395 = t254 * qJD(1);
t394 = t254 * qJD(3);
t257 = t268 * qJ(2);
t393 = t257 * qJD(1);
t392 = t268 * qJD(1);
t391 = t281 + qJD(5);
t387 = t501 / 0.2e1;
t385 = t92 * t538;
t383 = t359 * t431;
t382 = t290 * t539;
t381 = t530 * t431;
t380 = t558 * t430;
t378 = t359 * t282;
t377 = t530 * t420;
t376 = t530 * t419;
t134 = t359 * t538;
t375 = t253 * t395;
t374 = t294 * t539;
t372 = t530 * t429;
t370 = t547 / 0.2e1;
t368 = t530 * t507;
t365 = -t446 / 0.2e1;
t358 = t134 + t541;
t135 = t530 * t539;
t357 = -t135 - t541;
t355 = -qJD(1) * t274 - qJD(2);
t354 = qJD(6) + t538;
t279 = qJD(3) * t388;
t347 = 0.2e1 * t352;
t346 = t115 * qJD(4) + t283 * t359;
t345 = -t434 + t409;
t344 = t434 + t409;
t340 = t552 - t446;
t243 = -t273 * t388 - t277 * t501;
t303 = (t528 * t506 / 0.2e1 - t95 * t509) * pkin(3) + t95 * t514 + t528 * t513;
t325 = t555 / 0.2e1 - t564 / 0.2e1;
t5 = t303 + t325;
t339 = t5 * qJD(1) - t243 * qJD(3);
t336 = qJD(3) * t528 + qJD(4) * t42;
t334 = t462 + t551;
t330 = qJD(2) * t530 + qJD(4) * t43;
t328 = t354 * t294;
t327 = t354 * t297;
t326 = qJD(4) * t146 + t427;
t322 = -t554 / 0.2e1 - t547 / 0.2e1;
t321 = t387 - t361;
t320 = t271 * t515 + t546;
t319 = -t515 * t520 + t370;
t313 = t74 / 0.2e1 + t320;
t15 = t521 - t563 / 0.2e1 + t313 * t294;
t318 = -qJD(1) * t15 - t273 * t425;
t17 = t313 * t297;
t317 = -qJD(1) * t17 + t273 * t283;
t316 = t359 * (t423 + t427);
t312 = t79 / 0.2e1 + t319;
t10 = (-t302 + t529) * t294;
t310 = -t10 * qJD(1) - t279 * t297;
t12 = (t365 + t552 / 0.2e1 + t302) * t297;
t309 = -t12 * qJD(1) - t279 * t294;
t238 = t321 * t294;
t33 = t312 * t297;
t307 = qJ(5) * t282 - t33 * qJD(1) - t238 * qJD(3);
t239 = t321 * t297;
t32 = t312 * t294;
t306 = -qJ(5) * t421 - t32 * qJD(1) + t239 * qJD(3);
t289 = qJ(5) * qJD(5);
t269 = t288 * qJ(5);
t267 = t273 * qJD(5);
t262 = t270 * qJD(6);
t241 = t297 * t387 + (qJ(5) + t273) * t507;
t240 = 0.2e1 * t273 * t510;
t137 = t348 + t532;
t136 = t347 - t532;
t131 = t288 * t226;
t130 = t368 - t447 / 0.2e1;
t112 = t454 / 0.2e1 - t359 * t510;
t109 = t121 * qJD(6);
t108 = t130 * qJD(6);
t106 = t120 * qJD(6);
t104 = t350 * qJD(4);
t103 = t294 * t368 + t297 * t565;
t101 = -t107 - t419;
t54 = t535 * t530;
t49 = t370 + t554 / 0.2e1 + t322;
t47 = 0.2e1 * t360 * t530;
t45 = t502 / 0.2e1 - t546 + t530 * t513 - t322;
t31 = t300 + t323;
t22 = t297 * t319 + t508 * t79 - t486;
t21 = t563 - t495 / 0.2e1 + t319 * t294;
t18 = -t486 / 0.2e1 + t74 * t508 + t537 * t510 + t320 * t297;
t16 = t320 * t294 + t74 * t510 + 0.2e1 * t521;
t11 = t587 + (t302 + t529) * t297;
t9 = t586 + qJ(5) * t565 + (t365 - t302) * t294;
t6 = t303 - t325;
t13 = [0, 0, 0, 0, 0, t268 * qJD(2), t257 * qJD(2), t253 * t394, t172 * qJD(3), 0, 0, 0, t274 * t394, t274 * t252, -t54, -t288 * t562, 0, 0, 0, -qJD(3) * t80 + t422 * t530, -qJD(3) * t81 - t359 * t422, qJD(2) * t56, qJD(3) * t26 + qJD(4) * t28 + t398 * t530, qJD(3) * t27 + qJD(4) * t29 + qJD(5) * t558, qJD(2) * t25 + qJD(3) * t7 + qJD(4) * t8 - qJD(5) * t496, t290 * t316 + t373 * t559, -t132 * qJD(6) + 0.2e1 * t316 * t452, t288 * t573 + t359 * t376, t288 * t572 - t359 * t377, -t54, qJD(2) * t53 + qJD(3) * t1 + qJD(4) * t3 + qJD(6) * t20 + t284 * t558, -qJD(2) * t51 + qJD(3) * t2 + qJD(4) * t4 + qJD(6) * t19 + t285 * t558; 0, 0, 0, 0, 0, t392, t393, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, 0, t439, -t211, 0, qJD(3) * t45 + qJD(4) * t49 - qJD(5) * t228 + t478, 0, 0, 0, 0, 0, t104 + t108 + t469, -t471 + t580 + t582 + t583; 0, 0, 0, 0, 0, 0, 0, t375, t404, t252, -t394, 0, qJD(3) * t332 + t274 * t395, qJD(3) * t333 + t274 * t396, -t135, -t569, -t535, -t326, 0, -t336 - t466, -t560 - t465 (-t277 * t359 - t462) * qJD(3) + t31 * qJD(4) - t398, t336 + t477, t560 + t476, t489 + t45 * qJD(2) + (t273 * t95 + t277 * t528) * qJD(3) + t6 * qJD(4) + t42 * qJD(5), t103 * qJD(4) + t106 + (t283 * t297 + t382) * t530, t47 * qJD(4) + t173 + (-qJD(3) * t270 + t347) * t530, -t359 * t425 + t576 + t584, t108 + t346 + t585, t357, t498 + (-t297 * t334 + t587) * qJD(3) + t11 * qJD(4) - t123 * qJD(5) + t16 * qJD(6), t497 + t578 + (t294 * t334 + t586) * qJD(3) + t9 * qJD(4) + t112 * qJD(5) + t18 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t569, -t535, -qJD(3) * t146 - t423, 0, t381 - t389, -t560 - t383, t31 * qJD(3) + (-t552 + t505) * qJD(4) - t398, t389 + t475, t560 + t474, t485 + t49 * qJD(2) + t6 * qJD(3) + (-t555 + t564) * qJD(4) + t528 * qJD(5), t103 * qJD(3) + t106 + (t282 * t297 + t382) * t530, t47 * qJD(3) + t173 + (-qJD(4) * t270 + t347) * t530, qJD(3) * t574 - t359 * t421 + t584, qJD(3) * t115 + t378 + t585, -t358, t493 + t350 * qJD(2) + t11 * qJD(3) + (-t297 * t340 + t587) * qJD(4) + t574 * qJD(5) + t21 * qJD(6), t70 * t421 + t492 + t578 + t9 * qJD(3) + t22 * qJD(6) + (qJD(4) * t340 + t398) * t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t535, t134, t400, -t385 + t389, 0, 0, 0, 0, 0, -qJD(3) * t123 + t380 + t576, qJD(3) * t112 + t297 * t400 + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t523, -t561, qJD(3) * t571 + t354 * t548, t130 * qJD(3) - t328 * t359, -t131, qJD(2) * t130 + qJD(3) * t16 + qJD(4) * t21 + qJD(5) * t571 - qJD(6) * t24 + t479, qJD(2) * t571 + qJD(3) * t18 + qJD(4) * t22 + qJD(6) * t23 + t480; 0, 0, 0, 0, 0, -t392, -t393, 0, 0, 0, 0, 0, t394, t252, 0, 0, 0, 0, 0, t326, -t535, -t439, -t326, t535, -qJD(3) * t44 + qJD(4) * t48 - qJD(5) * t146 - t478, 0, 0, 0, 0, 0, -t109 + t346 - t469, qJD(3) * t548 + qJD(4) * t123 + t471 - t579; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t395, t396, 0, 0, 0, 0, 0, t538, -t539, 0, -t538, t539, -t442, 0, 0, 0, 0, 0, t374, t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t405, -t539, 0, -t405, t539, t441, 0, 0, 0, 0, 0, t413, t410; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t405, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t588; 0, 0, 0, 0, 0, 0, 0, -t375, -t404, 0, 0, 0, t355 * t254, t355 * t253, t135, t569, 0, -t211, 0, -t330 + t466, t542 + t465, qJD(4) * t30, t330 - t477, -t542 - t476, qJD(2) * t44 + qJD(4) * t5 + qJD(5) * t43 - t489, -qJD(4) * t102 - t382 * t530 + t106, t536, -qJD(6) * t110 + t582 - t584, t104 - t109 - t585, -t357, qJD(4) * t12 + qJD(6) * t15 - t294 * t542 - t498 - t581, -qJD(2) * t548 + qJD(4) * t10 + qJD(5) * t350 + qJD(6) * t17 - t497; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t395, -t396, 0, 0, 0, 0, 0, -t538, t539, 0, t538, -t539, t442, 0, 0, 0, 0, 0, -t374, -t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, -t281, 0, t280, t391, -qJD(4) * t243 + t267, -t373, t262, 0, 0, 0, t273 * t419 + t437, -t273 * t420 + t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t401, 0, -t526, -t281 - t279, t445, t526, t391 + t279 (-pkin(4) * t295 + qJ(5) * t506) * t499 + t267 + t339, -t373 - t433, t262, t409, t525, 0, t241 * qJD(6) - t309 + t437, t240 * qJD(6) - t310 + t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, t273 * t288 + t443, 0, 0, 0, 0, 0, -t345, t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t137, -t417 - t420, t101, t540, qJD(4) * t241 - t271 * t420 - t318, qJD(4) * t240 - t271 * t419 - t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t569, 0, t228 * qJD(3), 0, -t381 - t484, t542 + t383, -qJD(3) * t30, -t475 + t484, -t542 - t474, -qJD(2) * t48 - qJD(3) * t5 - t485, qJD(3) * t102 - t134 * t290 + t106, t536, -t377 - t584 - t583, -t376 - t585 - t524, t358, -qJD(2) * t115 - qJD(3) * t12 + qJD(6) * t32 - t493 + t581, -qJD(2) * t123 - qJD(3) * t10 + qJD(6) * t33 - t492; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t405, t539, 0, t405, -t539, -t441, 0, 0, 0, 0, 0, -t413, -t410; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t401, 0, t527, t279, -t445, -t527, qJD(5) - t279, t289 - t339, -t373 + t433, t262, -t409, -t525, 0, -t239 * qJD(6) + t284 + t309, t238 * qJD(6) + t285 + t310; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t289, -t373, t262, 0, 0, 0, qJ(5) * t419 + t284, -qJ(5) * t420 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, t269, 0, 0, 0, 0, 0, t344, t531; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t137, -t328, -t327, t540, t294 * t418 - t306, t297 * t418 - t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t400, t385 + t484, 0, 0, 0, 0, 0, -t380 + t579 - t582 + t583, -t524 + (-qJD(6) * t530 - t400) * t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t405, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, -qJ(5) * qJD(4) - qJD(3) * t273 - t443, 0, 0, 0, 0, 0, t345, -t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, -t269, 0, 0, 0, 0, 0, -t344, -t531; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t588, -t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t523, t561, t110 * qJD(3) + (-t359 * t429 + t282) * t530, t121 * qJD(3) + (t359 * t430 + t421) * t530, -t131, qJD(2) * t121 - qJD(3) * t15 - qJD(4) * t32 - qJD(5) * t575 - t479, qJD(2) * t575 - qJD(3) * t17 - qJD(4) * t33 + t285 * t530 - t480; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t136, t417, t107, -t540, qJD(4) * t239 + t318, -qJD(4) * t238 + t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t136, t530 * t430, t372, -t540, t306, t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t577, t372; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t13;