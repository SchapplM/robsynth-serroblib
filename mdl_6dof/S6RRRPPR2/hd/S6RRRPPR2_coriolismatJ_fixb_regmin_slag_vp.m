% Calculate minimal parameter regressor of coriolis matrix for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x30]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRRPPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:47
% EndTime: 2019-03-09 15:27:07
% DurationCPUTime: 11.40s
% Computational Cost: add. (10614->456), mult. (20557->583), div. (0->0), fcn. (23996->8), ass. (0->367)
t336 = qJD(2) + qJD(3);
t339 = sin(pkin(10));
t342 = sin(qJ(3));
t343 = sin(qJ(2));
t345 = cos(qJ(2));
t567 = cos(qJ(3));
t311 = -t342 * t345 - t343 * t567;
t586 = pkin(7) + pkin(8);
t317 = t586 * t343;
t318 = t586 * t345;
t371 = -t317 * t567 - t342 * t318;
t597 = t311 * qJ(4) + t371;
t614 = t339 * t597;
t309 = t342 * t343 - t345 * t567;
t312 = t567 * t318;
t494 = t342 * t317;
t596 = -t312 + t494;
t221 = t309 * qJ(4) + t596;
t340 = cos(pkin(10));
t629 = t340 * t221;
t640 = t629 / 0.2e1;
t202 = -t629 / 0.2e1;
t646 = t202 + t614 / 0.2e1;
t649 = -t614 / 0.2e1 + t640 + t646;
t659 = qJD(1) * t649;
t653 = 0.2e1 * t646;
t658 = qJD(3) * t653;
t657 = qJD(5) * t649;
t656 = t653 * qJD(5);
t344 = cos(qJ(6));
t568 = t344 / 0.2e1;
t499 = t340 * t309;
t505 = t339 * t311;
t269 = t499 - t505;
t636 = -t629 + t614;
t644 = -t269 * pkin(5) + t636;
t652 = t344 * t644;
t655 = t644 * t568 + t652 / 0.2e1;
t341 = sin(qJ(6));
t570 = -t341 / 0.2e1;
t654 = 0.2e1 * t644 * t570;
t598 = -t339 * t309 - t340 * t311;
t638 = t339 * t221 + t340 * t597;
t643 = -pkin(5) * t598 + t638;
t651 = t643 * t341;
t650 = t643 * t344;
t331 = -pkin(2) * t345 - pkin(1);
t282 = t309 * pkin(3) + t331;
t358 = -qJ(5) * t598 + t282;
t82 = (pkin(4) + pkin(9)) * t269 + t358;
t38 = t341 * t82 + t650;
t610 = t344 * t598;
t611 = t344 * t269;
t647 = t269 * t38 - t610 * t644 - t643 * t611;
t461 = qJD(6) * t341;
t613 = t341 * t598;
t622 = -t613 / 0.2e1;
t632 = 0.2e1 * t622;
t633 = t632 * qJD(1);
t645 = t633 - t461;
t642 = qJD(2) * t653;
t641 = t336 * t638;
t623 = -t611 / 0.2e1;
t631 = 0.2e1 * t623;
t639 = t336 * t631;
t496 = t341 * t269;
t621 = t613 / 0.2e1;
t624 = t622 + t621;
t635 = qJD(6) * t624;
t634 = qJD(6) * t632;
t337 = t341 ^ 2;
t338 = t344 ^ 2;
t320 = t337 - t338;
t619 = t269 ^ 2;
t158 = t320 * t619;
t602 = t336 * t344;
t612 = t341 * t602;
t364 = 0.2e1 * t612 * t269;
t626 = qJD(1) * t158 + t364;
t149 = (t338 / 0.2e1 - t337 / 0.2e1) * t269;
t110 = -qJD(1) * t149 + t612;
t608 = t598 * qJD(1);
t402 = qJD(6) + t608;
t618 = t598 ^ 2;
t584 = t598 / 0.2e1;
t561 = t339 * pkin(3);
t324 = qJ(5) + t561;
t616 = t324 * t269 / 0.2e1;
t498 = t340 * t342;
t442 = t567 * pkin(2);
t330 = t442 + pkin(3);
t504 = t339 * t330;
t296 = pkin(2) * t498 + t504;
t288 = qJ(5) + t296;
t581 = t288 / 0.2e1;
t615 = t269 * t581;
t603 = t336 * t269;
t325 = -t340 * pkin(3) - pkin(4);
t323 = -pkin(9) + t325;
t566 = pkin(2) * t342;
t322 = t339 * t566;
t303 = t340 * t442 - t322;
t433 = t567 * t339;
t302 = (t433 + t498) * pkin(2);
t577 = t302 / 0.2e1;
t373 = -t303 * t269 / 0.2e1 + t598 * t577;
t576 = -t324 / 0.2e1;
t350 = -t598 * (t576 + t581) + t373;
t295 = t330 * t340 - t322;
t289 = -pkin(4) - t295;
t287 = -pkin(9) + t289;
t583 = -t287 / 0.2e1;
t609 = t350 + t269 * (t323 / 0.2e1 + t583);
t463 = qJD(5) * t598;
t256 = t598 * qJD(4);
t607 = t336 * t371;
t605 = 0.2e1 * t341 * t611;
t406 = t598 * t336;
t405 = t336 * t341;
t599 = t288 + t324;
t388 = -t269 * t636 - t598 * t638;
t595 = qJD(1) * t388;
t594 = qJD(4) * t388;
t450 = t618 * qJD(1);
t387 = t619 + t618;
t591 = t387 * qJD(1);
t469 = qJD(1) * t344;
t429 = t341 * t469;
t590 = t149 * t336 + t429 * t619;
t162 = -0.2e1 * t269 * t429 + t320 * t336;
t589 = qJD(4) * t387;
t588 = pkin(4) / 0.2e1;
t585 = -qJ(5) / 0.2e1;
t582 = -t288 / 0.2e1;
t580 = t289 / 0.2e1;
t579 = t295 / 0.2e1;
t578 = -t296 / 0.2e1;
t399 = -t312 / 0.2e1;
t574 = -t325 / 0.2e1;
t335 = t343 * pkin(2);
t573 = t335 / 0.2e1;
t572 = t339 / 0.2e1;
t571 = -t340 / 0.2e1;
t569 = -t344 / 0.2e1;
t563 = t598 * pkin(4);
t562 = t311 * pkin(3);
t560 = pkin(3) * qJD(3);
t523 = t269 * qJ(5);
t125 = t523 - t562 + t563;
t114 = t125 + t335;
t262 = t598 * pkin(9);
t87 = t114 + t262;
t555 = t341 * t87;
t1 = (t652 - t555) * t598 + t647;
t559 = t1 * qJD(1);
t39 = t344 * t82 - t651;
t551 = t39 * t269;
t2 = t269 * t651 - t87 * t610 + t551;
t558 = t2 * qJD(1);
t93 = t125 + t262;
t554 = t341 * t93;
t3 = (t652 - t554) * t598 + t647;
t556 = t3 * qJD(1);
t4 = t496 * t643 - t93 * t610 + t551;
t550 = t4 * qJD(1);
t27 = t38 * t598 + t611 * t644;
t545 = qJD(1) * t27;
t28 = -t39 * t598 + t496 * t644;
t544 = qJD(1) * t28;
t113 = t269 * pkin(4) + t358;
t533 = t113 * t598;
t40 = -t114 * t269 - t533;
t543 = qJD(1) * t40;
t532 = t113 * t269;
t41 = -t114 * t598 + t532;
t542 = qJD(1) * t41;
t44 = -t125 * t269 - t533;
t539 = qJD(1) * t44;
t45 = -t125 * t598 + t532;
t538 = qJD(1) * t45;
t386 = -t619 + t618;
t63 = t386 * t341;
t537 = qJD(1) * t63;
t64 = t387 * t341;
t536 = qJD(1) * t64;
t65 = t386 * t344;
t535 = qJD(1) * t65;
t66 = t387 * t344;
t534 = qJD(1) * t66;
t13 = t113 * t114;
t531 = t13 * qJD(1);
t14 = t113 * t125;
t530 = t14 * qJD(1);
t522 = t269 * t323;
t521 = t288 * t598;
t520 = t289 * t269;
t29 = t282 * (t335 - t562);
t519 = t29 * qJD(1);
t518 = t295 * t269;
t517 = t296 * t598;
t30 = t282 * t562;
t516 = t30 * qJD(1);
t514 = t324 * t598;
t513 = t325 * t269;
t512 = t337 * t269;
t507 = t339 * t598;
t34 = -(t580 + t574) * t269 + t350;
t503 = t34 * qJD(1);
t500 = t340 * t269;
t347 = t518 / 0.2e1 - t517 / 0.2e1 + t373;
t361 = (-t507 / 0.2e1 + t500 / 0.2e1) * pkin(3);
t35 = t361 - t347;
t489 = t35 * qJD(1);
t301 = -t562 / 0.2e1;
t443 = t573 + t301;
t49 = -(t585 + t582) * t269 + (t588 - t289 / 0.2e1) * t598 + t443;
t488 = t49 * qJD(1);
t59 = t301 - (t585 + t576) * t269 + (t588 + t574) * t598;
t487 = t59 * qJD(1);
t376 = -t269 * t578 + t579 * t598;
t83 = t376 + t443;
t484 = t83 * qJD(1);
t418 = t269 * t570;
t419 = t496 / 0.2e1;
t143 = t419 - t418;
t153 = t623 + t611 / 0.2e1;
t483 = t153 * qJD(4) + t143 * qJD(5);
t482 = t631 * qJD(4);
t259 = t269 * qJD(4);
t481 = t341 * t259;
t298 = t303 * qJD(3);
t333 = qJD(5) * t341;
t478 = t341 * t298 + t333;
t334 = qJD(5) * t344;
t477 = t344 * t298 + t334;
t278 = t399 + t312 / 0.2e1;
t471 = qJD(1) * t278;
t470 = qJD(1) * t331;
t468 = qJD(1) * t345;
t467 = qJD(2) * t341;
t466 = qJD(2) * t344;
t465 = qJD(3) * t324;
t464 = qJD(3) * t331;
t462 = qJD(6) * t598;
t460 = qJD(6) * t344;
t362 = (-t269 * t572 + t571 * t598) * pkin(3);
t111 = t562 / 0.2e1 + t362;
t459 = t111 * qJD(1);
t141 = 0.2e1 * t621;
t128 = t141 * qJD(1);
t150 = 0.2e1 * t584 * t344;
t132 = t150 * qJD(1);
t456 = t631 * qJD(1);
t225 = t309 ^ 2 - t311 ^ 2;
t454 = t225 * qJD(1);
t252 = t309 * t335 - t311 * t331;
t453 = t252 * qJD(1);
t253 = -t309 * t331 - t311 * t335;
t452 = t253 * qJD(1);
t264 = -t499 / 0.2e1 + t505 / 0.2e1;
t451 = t264 * qJD(1);
t448 = t269 * qJD(1);
t321 = -t343 ^ 2 + t345 ^ 2;
t447 = t321 * qJD(1);
t446 = t343 * qJD(2);
t445 = t345 * qJD(2);
t444 = t298 + qJD(5);
t441 = pkin(1) * t343 * qJD(1);
t440 = pkin(1) * t468;
t432 = t113 * t608;
t431 = qJD(1) * t512;
t430 = t341 * t450;
t428 = t269 * t462;
t427 = t269 * t608;
t426 = t309 * t470;
t425 = t311 * t470;
t424 = t343 * t468;
t423 = t341 * t460;
t422 = t598 * t469;
t421 = -t521 / 0.2e1;
t420 = -t514 / 0.2e1;
t417 = t302 * t570;
t415 = t302 * t568;
t410 = t567 * qJD(2);
t409 = t567 * qJD(3);
t275 = t336 * t311;
t397 = t577 + t576 + t582;
t396 = t269 * t405;
t159 = t288 * t303 + t289 * t302;
t377 = -t638 * t577 + t636 * t303 / 0.2e1;
t349 = t580 * t636 + t581 * t638 + t377;
t378 = t574 * t636 + t576 * t638;
t6 = t349 + t378;
t393 = t6 * qJD(1) + t159 * qJD(2);
t160 = -t295 * t302 + t296 * t303;
t348 = t578 * t638 + t579 * t636 - t377;
t363 = (t571 * t636 + t572 * t638) * pkin(3);
t7 = t363 + t348;
t392 = t7 * qJD(1) - t160 * qJD(2);
t385 = t269 * t287 + t521;
t384 = t514 + t522;
t383 = qJD(2) * t288 - t659;
t52 = t202 + t640;
t382 = -qJD(1) * t52 - qJD(2) * t302;
t381 = qJD(2) * t303;
t9 = t341 * t609;
t379 = qJD(1) * t9 - t303 * t466;
t375 = -t522 / 0.2e1 + t420;
t374 = t583 * t598 + t615;
t372 = -t323 * t598 / 0.2e1 + t616;
t11 = t344 * t609;
t369 = -qJD(1) * t11 - t303 * t467;
t360 = t87 / 0.2e1 + t374;
t15 = t341 * t360;
t368 = -qJD(1) * t15 - t288 * t466;
t17 = t344 * t360;
t367 = -qJD(1) * t17 + t288 * t467;
t99 = qJD(6) * t264 - t448 * t598;
t366 = t301 + t563 / 0.2e1 + t523 / 0.2e1;
t359 = t93 / 0.2e1 + t372;
t279 = -qJ(5) + (t442 / 0.2e1 - pkin(3) / 0.2e1 - t330 / 0.2e1) * t339;
t357 = qJD(2) * t279 - t465 + t659;
t180 = t397 * t341;
t25 = t344 * t359;
t355 = -qJD(1) * t25 - qJD(2) * t180 + t341 * t465;
t181 = t397 * t344;
t23 = t341 * t359;
t354 = -qJD(1) * t23 + qJD(2) * t181 - t344 * t465;
t316 = t320 * qJD(6);
t297 = t302 * qJD(3);
t276 = t311 * t309 * qJD(1);
t274 = t336 * t309;
t273 = t561 / 0.2e1 + qJ(5) + t504 / 0.2e1 + (t498 + t433 / 0.2e1) * pkin(2);
t246 = t341 * t448;
t227 = 0.2e1 * t399 + t494;
t183 = t568 * t599 + t415;
t182 = t570 * t599 + t417;
t179 = t269 * qJD(5);
t156 = t336 * t264;
t140 = t150 * qJD(6);
t133 = t631 * qJD(5);
t131 = t149 * qJD(6);
t127 = -t132 - t460;
t126 = -t128 - t461;
t112 = t301 + t362;
t94 = t402 * t605;
t84 = -t376 + t443;
t62 = -t431 * t598 + t131;
t60 = t325 * t584 + t366 - t616;
t50 = t580 * t598 + t366 + t573 - t615;
t48 = -t140 - t535;
t47 = -qJD(6) * t141 - t537;
t46 = t131 + (t612 + t431) * t598;
t37 = (-qJD(6) + t608) * t605 - t320 * t406;
t36 = t361 + t347;
t33 = t421 - t520 / 0.2e1 + t420 - t513 / 0.2e1 + t373;
t32 = t396 + t535;
t31 = -t269 * t602 + t537 + t635;
t26 = t344 * t372 + t569 * t93 + t654;
t24 = -t554 / 0.2e1 + t372 * t341 + t655;
t18 = t344 * t374 + t569 * t87 + t654;
t16 = -t555 / 0.2e1 + t374 * t341 + t655;
t12 = t598 * t415 + t651 + (t421 + t375) * t344 + (t287 + t303) * t623;
t10 = -t287 * t418 + t288 * t621 + t303 * t419 - t375 * t341 + t598 * t417 + t650;
t8 = t363 - t348;
t5 = t349 - t378;
t19 = [0, 0, 0, t343 * t445, t321 * qJD(2), 0, 0, 0, -pkin(1) * t446, -pkin(1) * t445, t309 * t275, t336 * t225, 0, 0, 0, qJD(2) * t252 - t311 * t464, qJD(2) * t253 - t309 * t464, t589, qJD(2) * t29 - qJD(3) * t30 + t594, t589, qJD(2) * t40 + qJD(3) * t44 + t269 * t463, qJD(2) * t41 + qJD(3) * t45 + qJD(5) * t618, qJD(2) * t13 + qJD(3) * t14 - t113 * t463 + t594, t406 * t512 + t423 * t619, -qJD(6) * t158 + t364 * t598, t336 * t63 + t344 * t428, t336 * t65 - t341 * t428, -t598 * t603, qJD(2) * t1 + qJD(3) * t3 + qJD(4) * t66 + qJD(6) * t28 + t333 * t618, qJD(2) * t2 + qJD(3) * t4 - qJD(4) * t64 + qJD(6) * t27 + t334 * t618; 0, 0, 0, t424, t447, t445, -t446, 0, -pkin(7) * t445 - t441, pkin(7) * t446 - t440, t276, t454, -t274, t275, 0, qJD(2) * t596 + t227 * qJD(3) + t453, t452 - t607 (-t517 + t518) * qJD(2) + t36 * qJD(3), t519 + (-t295 * t636 + t296 * t638) * qJD(2) + t8 * qJD(3) + t84 * qJD(4) (-t520 - t521) * qJD(2) + t33 * qJD(3) - t179, qJD(2) * t636 + t543 + t658, t542 + t641, t531 + (t288 * t638 + t289 * t636) * qJD(2) + t5 * qJD(3) + t50 * qJD(4) + t656, t46, t37, t31, t32, t99, t559 + (-t344 * t385 + t651) * qJD(2) + t12 * qJD(3) + t133 + t16 * qJD(6), t558 + (t341 * t385 + t650) * qJD(2) + t10 * qJD(3) + t18 * qJD(6) + t483; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, t454, -t274, t275, 0, t227 * qJD(2) + qJD(3) * t596 - t425, -t426 - t607, t36 * qJD(2) + (t500 - t507) * t560, -t516 + t8 * qJD(2) + t112 * qJD(4) + (t339 * t638 - t340 * t636) * t560, t33 * qJD(2) + (-t513 - t514) * qJD(3) - t179, qJD(3) * t636 + t539 + t642, t538 + t641, t530 + t5 * qJD(2) + (t324 * t638 + t325 * t636) * qJD(3) + t60 * qJD(4) + t656, t46, t37, t31, t32, t99, t556 + t12 * qJD(2) + (-t344 * t384 + t651) * qJD(3) + t133 + t24 * qJD(6), t550 + t10 * qJD(2) + (t341 * t384 + t650) * qJD(3) + t26 * qJD(6) + t483; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t591, qJD(2) * t84 + qJD(3) * t112 + t595, t591, 0, 0, qJD(2) * t50 + qJD(3) * t60 + t595, 0, 0, 0, 0, 0, t534, t153 * t336 - t536 + t635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t603, t427, t450, -t432 + t642 + t658, 0, 0, 0, 0, 0, t430 + t635 + t639, t143 * t336 + t344 * t450; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t590, -t626, t336 * t624 + t402 * t611, -t402 * t496, t156, qJD(2) * t16 + qJD(3) * t24 + qJD(5) * t624 - qJD(6) * t39 + t544, qJD(2) * t18 + qJD(3) * t26 + qJD(4) * t624 + qJD(6) * t38 + t545; 0, 0, 0, -t424, -t447, 0, 0, 0, t441, t440, -t276, -t454, 0, 0, 0, qJD(3) * t278 - t453, -t452, -qJD(3) * t35, -qJD(3) * t7 - qJD(4) * t83 - t519, qJD(3) * t34, qJD(3) * t52 + t256 - t543, -t259 - t542, qJD(3) * t6 - qJD(4) * t49 - t531 - t657, t62, -t94, t47, t48, -t99, qJD(3) * t11 + qJD(6) * t15 - t481 - t559, -qJD(3) * t9 + qJD(6) * t17 + t482 - t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t566, -pkin(2) * t409, 0, t160 * qJD(3), 0, t297, t444, qJD(3) * t159 + qJD(5) * t288, -t423, t316, 0, 0, 0, t288 * t460 + t478, -t288 * t461 + t477; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t336 * t566 + t471 (-t410 - t409) * pkin(2), -t489 (-t302 * t340 + t303 * t339) * t560 - t392, t503, t297 - t382, t381 + t444 (t302 * t325 + t303 * t324) * qJD(3) + t273 * qJD(5) + t393, -t423, t316, 0, 0, 0, qJD(6) * t183 - t369 + t478, qJD(6) * t182 - t379 + t477; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t484, 0, t608, -t448, -t488, 0, 0, 0, 0, 0, -t246, t456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t336, qJD(3) * t273 + t383, 0, 0, 0, 0, 0, t405, t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, t162, t126, t127, -t451, qJD(3) * t183 - t287 * t461 - t368, qJD(3) * t182 - t287 * t460 - t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t276, -t454, 0, 0, 0, -qJD(2) * t278 + t425, t426, qJD(2) * t35, qJD(2) * t7 + qJD(4) * t111 + t516, -qJD(2) * t34, -qJD(2) * t52 + t256 - t539, -t259 - t538, -qJD(2) * t6 - qJD(4) * t59 - t530 - t657, t62, -t94, t47, t48, -t99, -qJD(2) * t11 + qJD(6) * t23 - t481 - t556, qJD(2) * t9 + qJD(6) * t25 + t482 - t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t566 - t471, pkin(2) * t410, t489, t392, -t503, t382, qJD(5) - t381, -qJD(5) * t279 - t393, -t423, t316, 0, 0, 0, -qJD(6) * t181 + t333 + t369, qJD(6) * t180 + t334 + t379; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t324 * qJD(5), -t423, t316, 0, 0, 0, t324 * t460 + t333, -t324 * t461 + t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, 0, t608, -t448, -t487, 0, 0, 0, 0, 0, -t246, t456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t336, -t357, 0, 0, 0, 0, 0, t405, t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, t162, t126, t127, -t451, -t323 * t461 - t354, -t323 * t460 - t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t591, qJD(2) * t83 - qJD(3) * t111 - t595, -t591, -t406, t603, qJD(2) * t49 + qJD(3) * t59 - t463 - t595, 0, 0, 0, 0, 0, -t140 + t396 - t534, t536 - t634 - t639; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t484, 0, -t608, t448, t488, 0, 0, 0, 0, 0, t246, -t456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t459, 0, -t608, t448, t487, 0, 0, 0, 0, 0, t246, -t456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t608, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t427, -t450, t336 * t649 + t256 + t432, 0, 0, 0, 0, 0, -t430 + t634 (-t450 - t462) * t344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t336, qJD(3) * t279 - t383, 0, 0, 0, 0, 0, -t405, -t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t336, t357, 0, 0, 0, 0, 0, -t405, -t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t608, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t645, -t402 * t344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t590, t626, t141 * t336 - t269 * t422, t150 * t336 + t341 * t427, t156, -qJD(2) * t15 - qJD(3) * t23 + qJD(4) * t150 - qJD(5) * t632 - t544, -qJD(2) * t17 - qJD(3) * t25 + qJD(4) * t632 + t334 * t598 - t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, -t162, t128, t132, t451, qJD(3) * t181 + t368, -qJD(3) * t180 + t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, -t162, t128, t132, t451, t354, t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t633; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t633, t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t19;
