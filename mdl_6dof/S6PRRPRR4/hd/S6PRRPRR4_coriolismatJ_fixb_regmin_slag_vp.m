% Calculate minimal parameter regressor of coriolis matrix for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x29]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRPRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:32
% EndTime: 2019-03-08 22:14:48
% DurationCPUTime: 7.66s
% Computational Cost: add. (3796->406), mult. (8899->607), div. (0->0), fcn. (9709->10), ass. (0->334)
t323 = sin(qJ(3));
t541 = pkin(8) - pkin(9);
t295 = t541 * t323;
t327 = cos(qJ(3));
t296 = t541 * t327;
t322 = sin(qJ(5));
t326 = cos(qJ(5));
t121 = t326 * t295 - t322 * t296;
t321 = sin(qJ(6));
t552 = t121 * t321;
t160 = -t552 / 0.2e1;
t161 = t552 / 0.2e1;
t563 = t160 + t161;
t325 = cos(qJ(6));
t551 = t121 * t325;
t162 = -t551 / 0.2e1;
t163 = t551 / 0.2e1;
t562 = t162 + t163;
t439 = qJD(3) - qJD(5);
t514 = qJ(4) * t323;
t529 = pkin(3) + pkin(4);
t259 = t327 * t529 + pkin(2) + t514;
t277 = t322 * t323 + t326 * t327;
t521 = pkin(5) * t277;
t309 = t322 * t327;
t485 = t326 * t323;
t535 = -t485 + t309;
t383 = pkin(10) * t535 + t521;
t332 = t259 + t383;
t361 = t322 * t295 + t326 * t296;
t545 = t325 * t361;
t46 = t321 * t332 + t545;
t561 = t535 * (t46 - t545);
t546 = t321 * t361;
t45 = -t325 * t332 + t546;
t560 = t535 * (t45 - t546);
t268 = t535 ^ 2;
t117 = -t277 ^ 2 + t268;
t553 = t117 * t325;
t559 = qJD(2) * t553;
t554 = t117 * t321;
t558 = qJD(2) * t554;
t557 = -pkin(10) / 0.2e1;
t550 = t117 * qJD(2);
t320 = sin(pkin(6));
t324 = sin(qJ(2));
t496 = t320 * t324;
t515 = cos(pkin(6));
t255 = t323 * t496 - t327 * t515;
t256 = t323 * t515 + t327 * t496;
t56 = t255 * t326 - t256 * t322;
t549 = t439 * t56;
t548 = t439 * t121;
t538 = t439 * t535;
t547 = t277 * t538;
t544 = t439 * t361;
t363 = t255 * t322 + t256 * t326;
t543 = t439 * t363;
t542 = -0.2e1 * t535;
t285 = qJ(4) * t322 + t326 * t529;
t281 = pkin(5) + t285;
t316 = t321 ^ 2;
t318 = t325 ^ 2;
t298 = t318 - t316;
t540 = t298 * t439;
t539 = t439 * t277;
t473 = qJD(2) * t535;
t184 = t277 * t473;
t260 = t485 / 0.2e1 - t309 / 0.2e1;
t482 = t260 * qJD(6) - t184;
t526 = -t281 / 0.2e1;
t406 = t285 / 0.2e1 + t526;
t286 = t326 * qJ(4) - t322 * t529;
t282 = -pkin(10) + t286;
t407 = -t282 / 0.2e1 + t286 / 0.2e1;
t537 = t277 * t406 - t535 * t407 - t521 / 0.2e1;
t405 = t316 / 0.2e1 - t318 / 0.2e1;
t139 = t405 * t535;
t454 = t139 * qJD(6);
t536 = t318 * t184 - t454;
t534 = -t326 * t277 - t322 * t535;
t531 = 0.2e1 * t281;
t469 = qJD(2) * t325;
t430 = t321 * t469;
t530 = -t139 * t439 - t268 * t430;
t528 = t277 / 0.2e1;
t527 = -t535 / 0.2e1;
t525 = -t321 / 0.2e1;
t524 = t321 / 0.2e1;
t523 = -t325 / 0.2e1;
t522 = t325 / 0.2e1;
t520 = pkin(5) * t535;
t519 = pkin(10) * t277;
t518 = t323 * pkin(3);
t328 = cos(qJ(2));
t495 = t320 * t328;
t507 = t363 * t321;
t77 = -t325 * t495 + t507;
t517 = t77 * t277;
t506 = t363 * t325;
t78 = t321 * t495 + t506;
t516 = t78 * t277;
t349 = t322 * t528 + t326 * t527;
t340 = t323 / 0.2e1 + t349;
t84 = t340 * t321;
t513 = qJD(2) * t84;
t87 = t340 * t325;
t512 = qJD(2) * t87;
t36 = t56 * t321;
t39 = t56 * t325;
t193 = t535 * t495;
t501 = t193 * t321;
t500 = t193 * t325;
t194 = t277 * t495;
t499 = t194 * t321;
t498 = t194 * t325;
t173 = t519 - t520;
t314 = t327 * qJ(4);
t265 = -t323 * t529 + t314;
t111 = -t173 + t265;
t494 = t321 * t111;
t140 = t321 * t277;
t142 = t321 * t535;
t489 = t325 * t111;
t145 = t325 * t277;
t146 = t325 * t535;
t58 = (t255 * t323 + t256 * t327 - t496) * t495;
t484 = t58 * qJD(1);
t414 = t495 / 0.2e1;
t389 = t327 * t414;
t391 = t323 * t414;
t480 = t322 * t391 + t326 * t389;
t415 = -t495 / 0.2e1;
t390 = t327 * t415;
t392 = t323 * t415;
t479 = t322 * t392 + t326 * t390;
t478 = t322 * t390 + t326 * t391;
t477 = t322 * t389 + t326 * t392;
t476 = qJD(2) * t139;
t151 = t298 * t268;
t475 = qJD(2) * t151;
t474 = qJD(2) * t277;
t472 = qJD(2) * t321;
t471 = qJD(2) * t323;
t470 = qJD(2) * t324;
t468 = qJD(2) * t327;
t467 = qJD(3) * qJ(4);
t466 = qJD(3) * t321;
t465 = qJD(3) * t325;
t464 = qJD(4) * t322;
t463 = qJD(4) * t326;
t462 = qJD(5) * t259;
t461 = qJD(5) * t321;
t460 = qJD(5) * t325;
t459 = qJD(6) * t321;
t458 = qJD(6) * t325;
t457 = qJD(6) * t326;
t453 = t140 * qJD(2);
t452 = t145 * qJD(2);
t377 = -pkin(3) * t327 - t514;
t287 = -pkin(2) + t377;
t292 = -t314 + t518;
t185 = t287 * t327 + t292 * t323;
t451 = t185 * qJD(2);
t186 = -t287 * t323 + t292 * t327;
t450 = t186 * qJD(2);
t449 = t260 * qJD(2);
t447 = t298 * qJD(6);
t317 = t323 ^ 2;
t319 = t327 ^ 2;
t299 = t319 - t317;
t446 = t299 * qJD(2);
t445 = t317 * qJD(2);
t444 = t322 * qJD(3);
t443 = t323 * qJD(3);
t442 = t323 * qJD(4);
t441 = t326 * qJD(3);
t313 = t327 * qJD(3);
t440 = t327 * qJD(4);
t437 = pkin(2) * t471;
t436 = pkin(2) * t468;
t435 = pkin(8) * t443;
t434 = pkin(8) * t313;
t432 = t325 * t496;
t431 = t318 * t473;
t429 = t277 * t442;
t428 = t277 * t459;
t427 = t277 * t458;
t426 = t287 * t471;
t425 = t320 * t470;
t424 = qJD(2) * t495;
t305 = t321 * t458;
t423 = t277 * t471;
t422 = t535 * t471;
t421 = t277 * t469;
t420 = t36 / 0.2e1;
t419 = t506 / 0.2e1;
t418 = -t499 / 0.2e1;
t417 = -t498 / 0.2e1;
t416 = t496 / 0.2e1;
t412 = -t146 / 0.2e1;
t404 = qJD(2) * (t317 + t319);
t402 = t439 * t321;
t400 = t439 * t325;
t291 = t439 * t326;
t399 = qJD(6) + t474;
t396 = t321 * t423;
t395 = t323 * t421;
t394 = t535 * t305;
t393 = t535 * t430;
t388 = t277 * t415;
t387 = t277 * t414;
t386 = t535 * t415;
t385 = t535 * t414;
t384 = -pkin(5) / 0.2e1 + t406;
t382 = 0.2e1 * t393;
t381 = -0.2e1 * t393;
t380 = t325 * t402;
t378 = t321 * t400;
t177 = -t500 / 0.2e1;
t339 = t363 * t527;
t330 = -t321 * t339 + t527 * t77;
t2 = t177 + t330;
t9 = t489 * t277 - t560;
t376 = t2 * qJD(1) + t9 * qJD(2);
t10 = -t494 * t277 - t561;
t174 = t501 / 0.2e1;
t329 = -t325 * t339 + t527 * t78;
t3 = t174 + t329;
t375 = t3 * qJD(1) + t10 * qJD(2);
t11 = t173 * t145 + t560;
t176 = t500 / 0.2e1;
t353 = (-t77 / 0.2e1 + t507 / 0.2e1) * t535;
t5 = t176 - t353;
t374 = t5 * qJD(1) + t11 * qJD(2);
t12 = -t173 * t140 + t561;
t175 = -t501 / 0.2e1;
t352 = (-t78 / 0.2e1 + t419) * t535;
t8 = t175 - t352;
t373 = t8 * qJD(1) + t12 * qJD(2);
t372 = t277 * t393;
t331 = -(t557 + t407) * t535 + (pkin(5) / 0.2e1 + t406) * t277;
t13 = t321 * t331;
t371 = qJD(2) * t13;
t16 = t325 * t331;
t370 = qJD(2) * t16;
t351 = t282 * t528 - t526 * t535;
t341 = t111 / 0.2e1 + t351;
t17 = t341 * t321 + t562;
t369 = qJD(2) * t17;
t19 = -t341 * t325 + t563;
t368 = qJD(2) * t19;
t350 = -t527 * t56 + t416;
t21 = t418 + t516 / 0.2e1 - t350 * t325;
t26 = t121 * t146 - t277 * t46;
t367 = qJD(1) * t21 - qJD(2) * t26;
t22 = t417 - t517 / 0.2e1 + t350 * t321;
t25 = -t121 * t142 + t277 * t45;
t366 = qJD(1) * t22 - qJD(2) * t25;
t75 = t259 * t535 + t265 * t277;
t92 = -t385 + t477;
t365 = qJD(1) * t92 - qJD(2) * t75;
t76 = t259 * t277 - t265 * t535;
t94 = t388 + t480;
t364 = qJD(1) * t94 - qJD(2) * t76;
t362 = t277 * t281 - t282 * t535;
t360 = t399 * t325;
t358 = qJD(3) * t281 - t463;
t357 = qJD(3) * t286 + t464;
t356 = qJD(5) * t286 + t464;
t355 = t519 / 0.2e1 - t520 / 0.2e1;
t354 = t314 / 0.2e1 - t518 / 0.2e1;
t95 = t387 + t479;
t348 = qJD(1) * t95 + t259 * t474;
t93 = -t386 + t478;
t347 = qJD(1) * t93 + t259 * t473;
t43 = (t292 / 0.2e1 + t354) * t495;
t346 = -t287 * t292 * qJD(2) + t43 * qJD(1);
t345 = t173 / 0.2e1 + t355;
t344 = t380 * t542;
t343 = -t255 * qJD(3) + t327 * t424;
t342 = t256 * qJD(3) + t323 * t424;
t304 = t321 * t444;
t338 = t322 * t461 - t325 * t457 - t304;
t307 = t325 * t444;
t337 = t321 * t457 + t322 * t460 - t307;
t133 = t384 * t321;
t29 = -t325 * t345 + t563;
t336 = pkin(5) * t461 - qJD(2) * t29 + qJD(3) * t133;
t134 = t384 * t325;
t27 = t321 * t345 + t562;
t335 = pkin(5) * t460 - qJD(2) * t27 + qJD(3) * t134;
t334 = qJD(3) * t377 + t440;
t306 = t323 * t468;
t290 = t439 * t322;
t200 = 0.2e1 * t394;
t199 = -0.2e1 * t394;
t198 = (-t324 * t468 - t328 * t443) * t320;
t197 = (-t313 * t328 + t323 * t470) * t320;
t150 = t439 * t260;
t136 = t531 * t522;
t135 = t531 * t524;
t128 = t321 * t145;
t120 = t382 - t540;
t119 = t381 + t540;
t99 = -t386 + t477;
t98 = -t385 + t478;
t97 = t388 + t479;
t96 = t387 + t480;
t86 = t323 * t522 - t325 * t349;
t85 = t321 * t349 + t323 * t525;
t80 = -t380 + t476;
t79 = t378 - t476;
t60 = t534 * t325;
t59 = t534 * t321;
t48 = t318 * t528 + (-t316 / 0.2e1 - t405) * t277;
t44 = t292 * t415 + t354 * t495;
t42 = -t363 * t523 + t419;
t41 = -t507 / 0.2e1 - t363 * t524;
t40 = t39 / 0.2e1 - t56 * t523;
t37 = -t525 * t56 + t420;
t30 = t173 * t522 - t325 * t355 + 0.2e1 * t160;
t28 = t173 * t525 + t321 * t355 + 0.2e1 * t162;
t24 = -t516 / 0.2e1 - t56 * t412 + t418 - t432 / 0.2e1;
t23 = t517 / 0.2e1 - t535 * t420 + t417 + t321 * t416;
t20 = 0.2e1 * t161 + t489 / 0.2e1 - t351 * t325;
t18 = 0.2e1 * t163 - t494 / 0.2e1 + t351 * t321;
t15 = pkin(10) * t412 + t537 * t325 - t546;
t14 = t142 * t557 + t537 * t321 + t545;
t7 = t174 - t352;
t6 = t177 - t353;
t4 = t175 + t329;
t1 = t176 + t330;
t31 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t425, -t424, 0, 0, 0, 0, 0, t198, t197, t198, t404 * t495, -t197, t484 + t44 * qJD(3) + (t287 * t470 + (pkin(8) * t404 + t442) * t328) * t320, 0, 0, 0, 0, 0, qJD(3) * t99 + qJD(5) * t98 - t277 * t425, qJD(3) * t96 + qJD(5) * t97 + t425 * t535, 0, 0, 0, 0, 0 ((-t432 - t499) * t277 - t193 * t142) * qJD(2) + t1 * qJD(3) + t6 * qJD(5) + t24 * qJD(6) (-(-t321 * t496 + t498) * t277 - t193 * t146) * qJD(2) + t4 * qJD(3) + t7 * qJD(5) + t23 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342, -t343, -t342, 0, t343, t44 * qJD(2) + (-pkin(3) * t256 - qJ(4) * t255) * qJD(3) + t256 * qJD(4), 0, 0, 0, 0, 0, qJD(2) * t99 - t543, qJD(2) * t96 - t549, 0, 0, 0, 0, 0, qJD(2) * t1 + qJD(5) * t42 + qJD(6) * t37 - t363 * t465, qJD(2) * t4 + qJD(5) * t41 + qJD(6) * t40 + t363 * t466; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t342, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t98 + t543, qJD(2) * t97 + t549, 0, 0, 0, 0, 0, qJD(2) * t6 + qJD(3) * t42 - qJD(6) * t36 - t363 * t460, qJD(2) * t7 + qJD(3) * t41 - qJD(6) * t39 + t363 * t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t24 + qJD(3) * t37 - qJD(5) * t36 - qJD(6) * t78, qJD(2) * t23 + qJD(3) * t40 - qJD(5) * t39 + qJD(6) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t43 - t484, 0, 0, 0, 0, 0, -qJD(3) * t92 - qJD(5) * t93, -qJD(3) * t94 - qJD(5) * t95, 0, 0, 0, 0, 0, qJD(3) * t2 + qJD(5) * t5 - qJD(6) * t21, qJD(3) * t3 + qJD(5) * t8 - qJD(6) * t22; 0, 0, 0, 0, t323 * t313, t299 * qJD(3), 0, 0, 0, -pkin(2) * t443, -pkin(2) * t313, -qJD(3) * t186 + t323 * t440, 0, -qJD(3) * t185 + qJD(4) * t317 (qJD(3) * t292 - t442) * t287, -t547, t439 * t117, 0, 0, 0, qJD(3) * t75 - t462 * t535 + t429, qJD(3) * t76 - t277 * t462 - t442 * t535, -t268 * t305 - t318 * t547, -qJD(6) * t151 - t277 * t344, t428 * t535 - t439 * t553, t427 * t535 + t439 * t554, t547, qJD(3) * t9 + qJD(5) * t11 + qJD(6) * t26 + t325 * t429, qJD(3) * t10 + qJD(5) * t12 + qJD(6) * t25 - t321 * t429; 0, 0, 0, 0, t306, t446, t313, -t443, 0, -t434 - t437, t435 - t436, -t434 - t450, t334, -t435 - t451, pkin(8) * t334 - t346, -t184, t550, -t539, t538, 0, -t365 - t544, -t364 - t548, qJD(5) * t128 - t454 + (-t321 * t465 - t431) * t277, t48 * qJD(5) + t199 + (-qJD(3) * t298 + t382) * t277, qJD(5) * t142 - t466 * t535 - t559, qJD(5) * t146 - t465 * t535 + t558, -t482 (t321 * t362 - t545) * qJD(3) + t59 * qJD(4) + t14 * qJD(5) + t20 * qJD(6) + t376 (t325 * t362 + t546) * qJD(3) + t60 * qJD(4) + t15 * qJD(5) + t18 * qJD(6) + t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t306, t313, t445, -t426 + t434, 0, 0, 0, 0, 0, t423, -t422, 0, 0, 0, 0, 0, qJD(3) * t59 + qJD(6) * t86 + t395, qJD(3) * t60 + qJD(6) * t85 - t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, -t550, t539, -t538, 0, -t347 + t544, -t348 + t548, qJD(3) * t128 + t454 + (-t321 * t460 + t431) * t277, t48 * qJD(3) + t200 + (-qJD(5) * t298 + t381) * t277, qJD(3) * t142 - t461 * t535 + t559, qJD(3) * t146 - t460 * t535 - t558, t482, t14 * qJD(3) + (t321 * t383 - t545) * qJD(5) + t30 * qJD(6) + t374, t15 * qJD(3) + (t325 * t383 + t546) * qJD(5) + t28 * qJD(6) + t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t530, t378 * t542 - t475, t399 * t142, t535 * t360, -t150, qJD(3) * t20 + qJD(4) * t86 + qJD(5) * t30 - qJD(6) * t46 - t367, qJD(3) * t18 + qJD(4) * t85 + qJD(5) * t28 + qJD(6) * t45 - t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t43, 0, 0, 0, 0, 0, qJD(2) * t92, qJD(2) * t94, 0, 0, 0, 0, 0, -qJD(2) * t2, -qJD(2) * t3; 0, 0, 0, 0, -t306, -t446, 0, 0, 0, t437, t436, t450, 0, t451, t346, t184, -t550, 0, 0, 0, t365, t364, t536, t199 - 0.2e1 * t372, -t427 + t559, t428 - t558, t482, qJD(5) * t13 + qJD(6) * t19 - t376, qJD(5) * t16 + qJD(6) * t17 - t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), 0, 0, 0, 0, 0, t356, -qJD(5) * t285 + t463, t305, t447, 0, 0, 0, -t281 * t459 + t325 * t356, -t281 * t458 - t321 * t356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t467, 0, 0, 0, 0, 0, t444, t441, 0, 0, 0, 0, 0, t307, -t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t439 * t286, -t285 * t439, -t305, -t447, 0, 0, 0, qJD(6) * t135 + t286 * t400 + t371, qJD(6) * t136 - t286 * t402 + t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t119, -t360, t399 * t321, t449, qJD(5) * t135 - t281 * t466 - t282 * t458 + t368, qJD(5) * t136 - t281 * t465 + t282 * t459 + t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t306, 0, -t445, t426, 0, 0, 0, 0, 0, -t423, t422, 0, 0, 0, 0, 0, -qJD(6) * t87 - t395, qJD(6) * t84 + t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t467, 0, 0, 0, 0, 0, -t290, -t291, 0, 0, 0, 0, 0, t337, -t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t291, 0, 0, 0, 0, 0, -t337, t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291 * t321 - t322 * t458 - t512, t291 * t325 + t322 * t459 + t513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t93, qJD(2) * t95, 0, 0, 0, 0, 0, -qJD(2) * t5, -qJD(2) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t550, 0, 0, 0, t347, t348, -t536, t200 + 0.2e1 * t372, qJD(6) * t145 - t559, -qJD(6) * t140 + t558, -t482, -qJD(3) * t13 + qJD(6) * t29 - t374, -qJD(3) * t16 + qJD(6) * t27 - t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t357, qJD(3) * t285 - t463, -t305, -t447, 0, 0, 0, -qJD(6) * t133 - t325 * t357 - t371, -qJD(6) * t134 + t321 * t357 - t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t444, -t441, 0, 0, 0, 0, 0, -t307, t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, t447, 0, 0, 0, -pkin(5) * t459, -pkin(5) * t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t120, t452 + t458, -t453 - t459, -t449, -pkin(10) * t458 - t336, pkin(10) * t459 - t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t21, qJD(2) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t530, -t344 + t475, -qJD(5) * t145 + (-t472 * t535 + t465) * t277, qJD(5) * t140 + (-t469 * t535 - t466) * t277, -t150, -qJD(3) * t19 + qJD(4) * t87 - qJD(5) * t29 + t367, -qJD(3) * t17 - qJD(4) * t84 - qJD(5) * t27 + t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t120, t421, -t277 * t472, -t449, qJD(5) * t133 + t321 * t358 - t368, qJD(5) * t134 + t325 * t358 - t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t321 * t441 + t512, -t325 * t441 - t513; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t119, -t452, t453, t449, t336, t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t31;
