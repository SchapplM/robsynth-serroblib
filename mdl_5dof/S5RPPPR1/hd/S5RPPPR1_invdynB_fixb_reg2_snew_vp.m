% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPPR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:51
% EndTime: 2022-01-20 09:13:03
% DurationCPUTime: 11.72s
% Computational Cost: add. (31526->552), mult. (75833->857), div. (0->0), fcn. (48937->10), ass. (0->397)
t564 = sin(qJ(1));
t566 = cos(qJ(1));
t533 = t566 * g(1) + t564 * g(2);
t567 = qJD(1) ^ 2;
t517 = -t567 * pkin(1) - t533;
t559 = sin(pkin(7));
t562 = cos(pkin(7));
t532 = t564 * g(1) - t566 * g(2);
t574 = qJDD(1) * pkin(1) + t532;
t456 = t559 * t517 - t562 * t574;
t457 = t562 * t517 + t559 * t574;
t588 = t559 * t456 + t562 * t457;
t401 = t562 * t456 - t559 * t457;
t615 = t566 * t401;
t662 = -t564 * t588 + t615;
t620 = t564 * t401;
t347 = t566 * t588 + t620;
t606 = t559 * qJDD(1);
t520 = t562 * t567 + t606;
t614 = g(3) - qJDD(2);
t492 = -qJ(2) * t520 + t562 * t614;
t604 = t562 * qJDD(1);
t521 = -t559 * t567 + t604;
t646 = -qJ(2) * t521 - t559 * t614;
t647 = t566 * t520 + t564 * t521;
t661 = -pkin(5) * t647 + t566 * t492 + t564 * t646;
t466 = -t564 * t520 + t566 * t521;
t660 = -pkin(5) * t466 - t564 * t492 + t566 * t646;
t556 = qJDD(1) * pkin(2);
t440 = -t567 * qJ(3) + qJDD(3) + t456 - t556;
t558 = sin(pkin(8));
t561 = cos(pkin(8));
t643 = pkin(3) * t561;
t581 = -qJ(4) * t558 - t643;
t613 = qJD(1) * t558;
t658 = -0.2e1 * qJD(4) * t613 + t581 * qJDD(1) + t440;
t557 = sin(pkin(9));
t560 = cos(pkin(9));
t563 = sin(qJ(5));
t565 = cos(qJ(5));
t577 = t557 * t565 + t560 * t563;
t573 = t577 * t558;
t487 = qJD(1) * t573;
t608 = qJDD(1) * t558;
t592 = t557 * t608;
t631 = t558 * t560;
t597 = t565 * t631;
t437 = -t487 * qJD(5) + qJDD(1) * t597 - t563 * t592;
t612 = t561 * qJD(1);
t538 = -qJD(5) + t612;
t640 = t487 * t538;
t657 = t437 + t640;
t441 = -t567 * pkin(2) + qJDD(1) * qJ(3) + t457;
t645 = 2 * qJD(3);
t656 = qJD(1) * t645 + t441;
t552 = t558 ^ 2;
t554 = t561 ^ 2;
t624 = t561 * t567;
t514 = (t552 + t554) * t624;
t590 = t561 * t604;
t472 = -t559 * t514 + t590;
t605 = t561 * qJDD(1);
t591 = t559 * t605;
t474 = t562 * t514 + t591;
t419 = t566 * t472 - t564 * t474;
t654 = t564 * t472 + t566 * t474;
t537 = -qJDD(5) + t605;
t489 = -t563 * t557 * t613 + qJD(1) * t597;
t639 = t489 * t487;
t572 = -t537 - t639;
t652 = t563 * t572;
t650 = t565 * t572;
t648 = (qJD(5) + t538) * t489;
t485 = t487 ^ 2;
t486 = t489 ^ 2;
t536 = t538 ^ 2;
t644 = pkin(3) * t558;
t638 = t538 * t563;
t637 = t538 * t565;
t636 = t552 * t567;
t426 = -t558 * t614 + t656 * t561;
t515 = t581 * qJD(1);
t395 = t515 * t612 + t426;
t575 = t658 * t560;
t576 = -pkin(4) * t561 - pkin(6) * t631;
t630 = t558 * t561;
t331 = t576 * qJDD(1) + (-t395 + (-pkin(4) * t552 * t560 + pkin(6) * t630) * t567) * t557 + t575;
t349 = t560 * t395 + t658 * t557;
t506 = t576 * qJD(1);
t550 = t557 ^ 2;
t535 = t550 * t636;
t332 = -pkin(4) * t535 - pkin(6) * t592 + t506 * t612 + t349;
t285 = -t565 * t331 + t563 * t332;
t286 = t563 * t331 + t565 * t332;
t247 = -t565 * t285 + t563 * t286;
t635 = t557 * t247;
t543 = t561 * t614;
t603 = qJDD(4) + t543;
t611 = t645 + t515;
t394 = (t611 * qJD(1) + t441) * t558 + t603;
t634 = t557 * t394;
t600 = t557 * t560 * t567;
t522 = t552 * t600;
t498 = -t522 + t605;
t633 = t557 * t498;
t499 = -t522 - t605;
t632 = t557 * t499;
t629 = t559 * t440;
t628 = t560 * t247;
t627 = t560 * t394;
t626 = t560 * t498;
t625 = t560 * t499;
t623 = t562 * t440;
t609 = qJDD(1) * t557;
t365 = -pkin(6) * t535 + (pkin(4) * t609 + t441 + (t506 * t560 + t611) * qJD(1)) * t558 + t603;
t622 = t563 * t365;
t429 = t537 - t639;
t621 = t563 * t429;
t617 = t565 * t365;
t616 = t565 * t429;
t607 = qJDD(1) * t560;
t553 = t560 ^ 2;
t601 = t553 * t636;
t599 = t557 * t624;
t598 = t558 * t639;
t596 = t558 * t624;
t595 = t560 * t624;
t594 = t561 * t639;
t593 = t557 * t607;
t589 = -t440 + t556;
t248 = t563 * t285 + t565 * t286;
t425 = t656 * t558 + t543;
t364 = t558 * t425 + t561 * t426;
t477 = -t564 * t532 - t566 * t533;
t551 = t558 * t552;
t585 = t551 * t600;
t584 = t557 * t595;
t583 = t558 * t590;
t529 = t566 * qJDD(1) - t564 * t567;
t582 = -pkin(5) * t529 - t564 * g(3);
t580 = t552 * t584;
t348 = t557 * t395 - t575;
t300 = -t560 * t348 + t557 * t349;
t301 = t557 * t348 + t560 * t349;
t363 = t561 * t425 - t558 * t426;
t481 = t520 * t630;
t482 = -t559 * t596 + t583;
t579 = t566 * t481 + t564 * t482;
t578 = t564 * t481 - t566 * t482;
t476 = t566 * t532 - t564 * t533;
t571 = qJDD(1) * t573;
t570 = t577 * t608;
t547 = t554 * t567;
t545 = t554 * qJDD(1);
t544 = t552 * qJDD(1);
t528 = t564 * qJDD(1) + t566 * t567;
t525 = t547 - t636;
t524 = t547 + t636;
t519 = t545 - t544;
t518 = t545 + t544;
t513 = -t547 - t601;
t512 = t547 - t601;
t511 = (t554 * t558 + t551) * t567;
t510 = -t535 - t547;
t509 = t535 - t547;
t502 = -pkin(5) * t528 + t566 * g(3);
t501 = t535 - t601;
t500 = t535 + t601;
t497 = (t599 - t607) * t558;
t496 = (t599 + t607) * t558;
t495 = (t595 - t609) * t558;
t494 = (t595 + t609) * t558;
t484 = (-t550 - t553) * t596;
t483 = (qJDD(1) * t553 + t584) * t558;
t480 = (t553 * t624 - t593) * t558;
t479 = (t550 * t624 + t593) * t558;
t478 = (qJDD(1) * t550 - t584) * t558;
t473 = t562 * t511 + t558 * t606;
t470 = t559 * t511 - t558 * t604;
t463 = t562 * t519 - t559 * t525;
t462 = t562 * t518 - t559 * t524;
t461 = t559 * t519 + t562 * t525;
t460 = t559 * t518 + t562 * t524;
t459 = -t486 + t536;
t458 = t485 - t536;
t454 = t561 * t483 + t585;
t453 = t561 * t478 - t585;
t452 = -t559 * t484 - t583;
t451 = t562 * t484 - t558 * t591;
t450 = -t486 - t536;
t449 = -t557 * t512 + t625;
t448 = -t557 * t513 + t626;
t447 = t560 * t510 - t632;
t446 = t560 * t509 + t633;
t445 = -t560 * t512 - t632;
t444 = t560 * t513 + t633;
t443 = t557 * t510 + t625;
t442 = -t557 * t509 + t626;
t438 = -t486 + t485;
t436 = -t489 * qJD(5) - t571;
t435 = -t560 * t494 - t557 * t497;
t434 = t560 * t495 - t557 * t496;
t433 = -t557 * t494 + t560 * t497;
t432 = -t557 * t495 - t560 * t496;
t427 = -t536 - t485;
t420 = -t564 * t470 + t566 * t473;
t418 = t566 * t470 + t564 * t473;
t417 = (t487 * t565 - t489 * t563) * t538;
t416 = (t487 * t563 + t489 * t565) * t538;
t415 = t561 * t449 - t558 * t497;
t414 = t561 * t448 + t558 * t496;
t413 = t561 * t447 - t558 * t495;
t412 = t561 * t446 - t558 * t494;
t411 = t558 * t448 - t561 * t496;
t410 = t558 * t447 + t561 * t495;
t409 = -t564 * t460 + t566 * t462;
t408 = t566 * t460 + t564 * t462;
t407 = -t485 - t486;
t406 = t562 * t454 - t559 * t480;
t405 = t562 * t453 - t559 * t479;
t404 = t559 * t454 + t562 * t480;
t403 = t559 * t453 + t562 * t479;
t400 = t561 * t435 - t558 * t500;
t399 = t561 * t434 - t558 * t501;
t396 = t558 * t435 + t561 * t500;
t393 = -t437 + t640;
t390 = -t571 - t648;
t389 = t570 + t648;
t388 = (qJD(5) - t538) * t489 + t570;
t386 = pkin(1) * t614 + qJ(2) * t588;
t385 = t565 * t437 + t489 * t638;
t384 = t563 * t437 - t489 * t637;
t383 = -t563 * t436 - t487 * t637;
t382 = t565 * t436 - t487 * t638;
t381 = t565 * t458 + t621;
t380 = -t563 * t459 + t650;
t379 = t563 * t458 - t616;
t378 = t565 * t459 + t652;
t377 = -t563 * t450 + t616;
t376 = t565 * t450 + t621;
t375 = t565 * t427 - t652;
t374 = t563 * t427 + t650;
t373 = t562 * t415 - t559 * t445;
t372 = t562 * t414 + t559 * t444;
t371 = t562 * t413 + t559 * t443;
t370 = t562 * t412 - t559 * t442;
t369 = t559 * t415 + t562 * t445;
t368 = t559 * t414 - t562 * t444;
t367 = t559 * t413 - t562 * t443;
t366 = t559 * t412 + t562 * t442;
t361 = -qJ(4) * t444 + t627;
t360 = -qJ(4) * t443 + t634;
t359 = t562 * t400 + t559 * t433;
t358 = t562 * t399 - t559 * t432;
t357 = t559 * t400 - t562 * t433;
t356 = t559 * t399 + t562 * t432;
t355 = -qJ(2) * t470 - t559 * t426 + t561 * t623;
t354 = -qJ(2) * t472 - t559 * t425 + t558 * t623;
t353 = qJ(2) * t473 + t562 * t426 + t561 * t629;
t352 = -qJ(2) * t474 + t562 * t425 + t558 * t629;
t351 = -t557 * t416 + t560 * t417;
t350 = -t560 * t416 - t557 * t417;
t345 = t561 * t351 - t558 * t537;
t344 = -qJ(2) * t460 + t562 * t363;
t343 = qJ(2) * t462 + t559 * t363;
t342 = t565 * t390 - t563 * t393;
t341 = -t565 * t388 - t563 * t657;
t340 = t563 * t390 + t565 * t393;
t339 = -t563 * t388 + t565 * t657;
t338 = t562 * t364 + t629;
t337 = t559 * t364 - t623;
t336 = -t557 * t384 + t560 * t385;
t335 = -t557 * t382 + t560 * t383;
t334 = -t560 * t384 - t557 * t385;
t333 = -t560 * t382 - t557 * t383;
t330 = -pkin(3) * t444 + t349;
t329 = -pkin(3) * t443 + t348;
t325 = -t557 * t379 + t560 * t381;
t324 = -t557 * t378 + t560 * t380;
t323 = -t560 * t379 - t557 * t381;
t322 = -t560 * t378 - t557 * t380;
t321 = -t557 * t376 + t560 * t377;
t320 = t560 * t376 + t557 * t377;
t319 = -pkin(2) * t410 - pkin(3) * t495 - qJ(4) * t447 + t627;
t318 = -pkin(2) * t411 + pkin(3) * t496 - qJ(4) * t448 - t634;
t317 = -pkin(6) * t376 + t617;
t316 = t561 * t336 + t598;
t315 = t561 * t335 - t598;
t314 = -t557 * t374 + t560 * t375;
t313 = t560 * t374 + t557 * t375;
t312 = -t564 * t368 + t566 * t372;
t311 = -t564 * t367 + t566 * t371;
t310 = t566 * t368 + t564 * t372;
t309 = t566 * t367 + t564 * t371;
t308 = -pkin(6) * t374 + t622;
t307 = -t564 * t357 + t566 * t359;
t306 = t566 * t357 + t564 * t359;
t305 = t561 * t325 - t558 * t389;
t304 = t561 * t324 - t558 * t393;
t303 = t561 * t321 + t558 * t657;
t302 = t558 * t321 - t561 * t657;
t299 = -pkin(4) * t657 + pkin(6) * t377 + t622;
t298 = t562 * t345 - t559 * t350;
t297 = t559 * t345 + t562 * t350;
t296 = -pkin(4) * t388 + pkin(6) * t375 - t617;
t295 = t561 * t314 + t558 * t388;
t294 = t558 * t314 - t561 * t388;
t293 = -t557 * t340 + t560 * t342;
t292 = -t557 * t339 + t560 * t341;
t291 = t560 * t340 + t557 * t342;
t290 = -t560 * t339 - t557 * t341;
t289 = -t564 * t337 + t566 * t338;
t288 = t566 * t337 + t564 * t338;
t287 = -qJ(4) * t433 - t300;
t283 = t561 * t301 + t558 * t394;
t282 = t558 * t301 - t561 * t394;
t281 = t561 * t292 - t558 * t438;
t280 = -qJ(3) * t411 - t558 * t330 + t561 * t361;
t279 = -qJ(3) * t410 - t558 * t329 + t561 * t360;
t278 = t561 * t293 + t558 * t407;
t277 = t558 * t293 - t561 * t407;
t276 = t562 * t316 - t559 * t334;
t275 = t562 * t315 - t559 * t333;
t274 = t559 * t316 + t562 * t334;
t273 = t559 * t315 + t562 * t333;
t272 = -qJ(2) * t337 - (pkin(2) * t559 - qJ(3) * t562) * t363;
t271 = -pkin(2) * t396 - pkin(3) * t500 - qJ(4) * t435 - t301;
t270 = t562 * t305 - t559 * t323;
t269 = t562 * t304 - t559 * t322;
t268 = t559 * t305 + t562 * t323;
t267 = t559 * t304 + t562 * t322;
t266 = -qJ(3) * t396 + t561 * t287 + t433 * t644;
t265 = -pkin(3) * t291 - pkin(4) * t340;
t264 = t562 * t303 + t559 * t320;
t263 = t559 * t303 - t562 * t320;
t262 = qJ(2) * t338 - (-pkin(2) * t562 - qJ(3) * t559 - pkin(1)) * t363;
t261 = t562 * t295 + t559 * t313;
t260 = t559 * t295 - t562 * t313;
t259 = -pkin(3) * t320 - pkin(4) * t376 + t286;
t258 = t562 * t283 + t559 * t300;
t257 = t559 * t283 - t562 * t300;
t256 = -qJ(2) * t368 + t562 * t280 - t559 * t318;
t255 = -qJ(2) * t367 + t562 * t279 - t559 * t319;
t254 = t562 * t281 - t559 * t290;
t253 = t559 * t281 + t562 * t290;
t252 = -pkin(3) * t313 - pkin(4) * t374 + t285;
t251 = -qJ(4) * t320 - t557 * t299 + t560 * t317;
t250 = t562 * t278 + t559 * t291;
t249 = t559 * t278 - t562 * t291;
t246 = -pkin(2) * t282 + pkin(3) * t394 - qJ(4) * t301;
t245 = -pkin(1) * t411 + qJ(2) * t372 + t559 * t280 + t562 * t318;
t244 = -pkin(1) * t410 + qJ(2) * t371 + t559 * t279 + t562 * t319;
t243 = -qJ(4) * t313 - t557 * t296 + t560 * t308;
t242 = -pkin(4) * t365 + pkin(6) * t248;
t241 = -pkin(6) * t340 - t247;
t240 = -pkin(4) * t407 + pkin(6) * t342 + t248;
t239 = -qJ(3) * t282 + (-qJ(4) * t561 + t644) * t300;
t238 = -t564 * t263 + t566 * t264;
t237 = t566 * t263 + t564 * t264;
t236 = -qJ(2) * t357 + t562 * t266 - t559 * t271;
t235 = -pkin(2) * t302 + pkin(3) * t657 - qJ(4) * t321 - t560 * t299 - t557 * t317;
t234 = -t564 * t260 + t566 * t261;
t233 = t566 * t260 + t564 * t261;
t232 = -pkin(1) * t396 + qJ(2) * t359 + t559 * t266 + t562 * t271;
t231 = -pkin(2) * t294 + pkin(3) * t388 - qJ(4) * t314 - t560 * t296 - t557 * t308;
t230 = -t564 * t257 + t566 * t258;
t229 = t566 * t257 + t564 * t258;
t228 = -t564 * t249 + t566 * t250;
t227 = t566 * t249 + t564 * t250;
t226 = t560 * t248 - t635;
t225 = t557 * t248 + t628;
t224 = t561 * t226 + t558 * t365;
t223 = t558 * t226 - t561 * t365;
t222 = -qJ(3) * t302 + t561 * t251 - t558 * t259;
t221 = -qJ(3) * t294 + t561 * t243 - t558 * t252;
t220 = -qJ(4) * t291 - t557 * t240 + t560 * t241;
t219 = -pkin(3) * t225 - pkin(4) * t247;
t218 = -qJ(2) * t257 + t562 * t239 - t559 * t246;
t217 = -pkin(2) * t277 + pkin(3) * t407 - qJ(4) * t293 - t560 * t240 - t557 * t241;
t216 = -pkin(1) * t282 + qJ(2) * t258 + t559 * t239 + t562 * t246;
t215 = -qJ(3) * t277 + t561 * t220 - t558 * t265;
t214 = -pkin(6) * t628 - qJ(4) * t225 - t557 * t242;
t213 = t562 * t224 + t559 * t225;
t212 = t559 * t224 - t562 * t225;
t211 = -qJ(2) * t263 + t562 * t222 - t559 * t235;
t210 = -qJ(2) * t260 + t562 * t221 - t559 * t231;
t209 = -pkin(1) * t302 + qJ(2) * t264 + t559 * t222 + t562 * t235;
t208 = -pkin(1) * t294 + qJ(2) * t261 + t559 * t221 + t562 * t231;
t207 = -pkin(2) * t223 + pkin(3) * t365 + pkin(6) * t635 - qJ(4) * t226 - t560 * t242;
t206 = -t564 * t212 + t566 * t213;
t205 = t566 * t212 + t564 * t213;
t204 = -qJ(2) * t249 + t562 * t215 - t559 * t217;
t203 = -pkin(1) * t277 + qJ(2) * t250 + t559 * t215 + t562 * t217;
t202 = -qJ(3) * t223 + t561 * t214 - t558 * t219;
t201 = -qJ(2) * t212 + t562 * t202 - t559 * t207;
t200 = -pkin(1) * t223 + qJ(2) * t213 + t559 * t202 + t562 * t207;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t528, -t529, 0, t477, 0, 0, 0, 0, 0, 0, -t647, -t466, 0, t347, 0, 0, 0, 0, 0, 0, -t654, t420, t409, t289, 0, 0, 0, 0, 0, 0, t311, t312, t307, t230, 0, 0, 0, 0, 0, 0, t234, t238, t228, t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t529, -t528, 0, t476, 0, 0, 0, 0, 0, 0, t466, -t647, 0, -t662, 0, 0, 0, 0, 0, 0, t419, t418, t408, t288, 0, 0, 0, 0, 0, 0, t309, t310, t306, t229, 0, 0, 0, 0, 0, 0, t233, t237, t227, t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t614, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t363, 0, 0, 0, 0, 0, 0, t410, t411, t396, t282, 0, 0, 0, 0, 0, 0, t294, t302, t277, t223; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t529, 0, -t528, 0, t582, -t502, -t476, -pkin(5) * t476, 0, 0, t466, 0, -t647, 0, t660, -t661, t662, pkin(5) * t662 + qJ(2) * t615 - t564 * t386, -t578, -t564 * t461 + t566 * t463, t420, t578, t654, 0, -pkin(5) * t419 - t564 * t352 + t566 * t354, -pkin(5) * t418 - t564 * t353 + t566 * t355, -pkin(5) * t408 - t564 * t343 + t566 * t344, -pkin(5) * t288 - t564 * t262 + t566 * t272, -t564 * t404 + t566 * t406, -t564 * t356 + t566 * t358, -t564 * t369 + t566 * t373, -t564 * t403 + t566 * t405, -t564 * t366 + t566 * t370, -t564 * t451 + t566 * t452, -pkin(5) * t309 - t564 * t244 + t566 * t255, -pkin(5) * t310 - t564 * t245 + t566 * t256, -pkin(5) * t306 - t564 * t232 + t566 * t236, -pkin(5) * t229 - t564 * t216 + t566 * t218, -t564 * t274 + t566 * t276, -t564 * t253 + t566 * t254, -t564 * t267 + t566 * t269, -t564 * t273 + t566 * t275, -t564 * t268 + t566 * t270, -t564 * t297 + t566 * t298, -pkin(5) * t233 - t564 * t208 + t566 * t210, -pkin(5) * t237 - t564 * t209 + t566 * t211, -pkin(5) * t227 - t564 * t203 + t566 * t204, -pkin(5) * t205 - t564 * t200 + t566 * t201; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t528, 0, t529, 0, t502, t582, t477, pkin(5) * t477, 0, 0, t647, 0, t466, 0, t661, t660, t347, pkin(5) * t347 + qJ(2) * t620 + t566 * t386, t579, t566 * t461 + t564 * t463, t418, -t579, -t419, 0, -pkin(5) * t654 + t566 * t352 + t564 * t354, pkin(5) * t420 + t566 * t353 + t564 * t355, pkin(5) * t409 + t566 * t343 + t564 * t344, pkin(5) * t289 + t566 * t262 + t564 * t272, t566 * t404 + t564 * t406, t566 * t356 + t564 * t358, t566 * t369 + t564 * t373, t566 * t403 + t564 * t405, t566 * t366 + t564 * t370, t566 * t451 + t564 * t452, pkin(5) * t311 + t566 * t244 + t564 * t255, pkin(5) * t312 + t566 * t245 + t564 * t256, pkin(5) * t307 + t566 * t232 + t564 * t236, pkin(5) * t230 + t566 * t216 + t564 * t218, t566 * t274 + t564 * t276, t566 * t253 + t564 * t254, t566 * t267 + t564 * t269, t566 * t273 + t564 * t275, t566 * t268 + t564 * t270, t566 * t297 + t564 * t298, pkin(5) * t234 + t566 * t208 + t564 * t210, pkin(5) * t238 + t566 * t209 + t564 * t211, pkin(5) * t228 + t566 * t203 + t564 * t204, pkin(5) * t206 + t566 * t200 + t564 * t201; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t532, t533, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t521 - t456, -pkin(1) * t520 - t457, 0, -pkin(1) * t401, t544, 0.2e1 * t558 * t605, 0, t545, 0, 0, pkin(1) * t472 - qJ(3) * t514 + t561 * t589, pkin(1) * t470 + qJ(3) * t511 - t558 * t589, pkin(1) * t460 + pkin(2) * t524 + qJ(3) * t518 + t364, pkin(1) * t337 - pkin(2) * t440 + qJ(3) * t364, t558 * t483 - t580, t558 * t434 + t561 * t501, t558 * t449 + t561 * t497, t558 * t478 + t580, t558 * t446 + t561 * t494, t545, pkin(1) * t367 - pkin(2) * t443 + qJ(3) * t413 + t561 * t329 + t558 * t360, pkin(1) * t368 - pkin(2) * t444 + qJ(3) * t414 + t561 * t330 + t558 * t361, pkin(1) * t357 + qJ(3) * t400 + t558 * t287 + (-pkin(2) - t643) * t433, pkin(1) * t257 + qJ(3) * t283 + (-pkin(2) + t581) * t300, t558 * t336 - t594, t558 * t292 + t561 * t438, t558 * t324 + t561 * t393, t558 * t335 + t594, t558 * t325 + t561 * t389, t558 * t351 + t561 * t537, pkin(1) * t260 - pkin(2) * t313 + qJ(3) * t295 + t558 * t243 + t561 * t252, pkin(1) * t263 - pkin(2) * t320 + qJ(3) * t303 + t558 * t251 + t561 * t259, pkin(1) * t249 - pkin(2) * t291 + qJ(3) * t278 + t558 * t220 + t561 * t265, pkin(1) * t212 - pkin(2) * t225 + qJ(3) * t224 + t558 * t214 + t561 * t219;];
tauB_reg = t1;
