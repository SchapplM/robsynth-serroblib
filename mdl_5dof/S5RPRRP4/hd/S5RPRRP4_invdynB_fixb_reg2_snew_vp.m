% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:34
% EndTime: 2020-01-03 11:50:46
% DurationCPUTime: 10.72s
% Computational Cost: add. (29115->512), mult. (71637->741), div. (0->0), fcn. (49311->8), ass. (0->388)
t686 = 2 * qJD(2);
t590 = sin(qJ(1));
t593 = cos(qJ(1));
t559 = t590 * g(2) - t593 * g(3);
t594 = qJD(1) ^ 2;
t543 = -t594 * pkin(1) + qJDD(1) * qJ(2) - t559;
t692 = (qJD(1) * t686) + t543;
t588 = sin(qJ(4));
t589 = sin(qJ(3));
t591 = cos(qJ(4));
t592 = cos(qJ(3));
t586 = sin(pkin(8));
t636 = qJD(1) * t586;
t526 = (-t588 * t592 - t589 * t591) * t636;
t619 = t589 * t636;
t535 = t592 * t586 * qJDD(1) - qJD(3) * t619;
t630 = qJDD(1) * t589;
t635 = qJD(1) * t592;
t602 = qJD(3) * t635 + t630;
t597 = t602 * t586;
t468 = t526 * qJD(4) + t591 * t535 - t588 * t597;
t587 = cos(pkin(8));
t633 = t587 * qJD(1);
t569 = -qJD(3) + t633;
t562 = -qJD(4) + t569;
t666 = t526 * t562;
t693 = t468 - t666;
t629 = t587 * qJDD(1);
t568 = -qJDD(3) + t629;
t558 = -qJDD(4) + t568;
t618 = t586 * t635;
t528 = -t588 * t619 + t591 * t618;
t664 = t528 * t526;
t598 = -t558 + t664;
t691 = t598 * pkin(4);
t655 = t588 * t598;
t646 = t591 * t598;
t609 = t569 * t619;
t500 = t535 + t609;
t490 = t592 * t500;
t525 = t528 ^ 2;
t557 = t562 ^ 2;
t488 = -t525 - t557;
t473 = t558 + t664;
t656 = t588 * t473;
t429 = t591 * t488 + t656;
t647 = t591 * t473;
t430 = -t588 * t488 + t647;
t386 = t592 * t429 + t589 * t430;
t690 = -pkin(2) * t386 - pkin(3) * t429;
t524 = t526 ^ 2;
t472 = -t557 - t524;
t420 = t588 * t472 + t646;
t421 = t591 * t472 - t655;
t376 = t592 * t420 + t589 * t421;
t689 = -pkin(2) * t376 - pkin(3) * t420;
t612 = t588 * t535 + t591 * t597;
t432 = (qJD(4) + t562) * t528 + t612;
t467 = -t528 * qJD(4) - t612;
t494 = -t562 * pkin(4) - t528 * qJ(5);
t688 = -t467 * pkin(4) - t524 * qJ(5) + t528 * t494 + qJDD(5);
t566 = t569 ^ 2;
t687 = t589 ^ 2;
t683 = pkin(2) * t586;
t682 = pkin(2) * t587;
t436 = -t666 - t468;
t391 = -t432 * t588 + t591 * t436;
t393 = -t432 * t591 - t588 * t436;
t344 = -t589 * t391 + t592 * t393;
t469 = -t524 - t525;
t331 = t587 * t344 + t586 * t469;
t342 = t592 * t391 + t589 * t393;
t291 = -t593 * t331 - t590 * t342;
t679 = pkin(5) * t291;
t377 = -t589 * t420 + t592 * t421;
t431 = (qJD(4) - t562) * t528 + t612;
t349 = t587 * t377 + t586 * t431;
t315 = -t593 * t349 - t590 * t376;
t678 = pkin(5) * t315;
t388 = -t589 * t429 + t592 * t430;
t355 = t587 * t388 + t586 * t693;
t320 = -t593 * t355 - t590 * t386;
t677 = pkin(5) * t320;
t676 = pkin(6) * t342;
t675 = pkin(6) * t376;
t674 = pkin(6) * t386;
t673 = pkin(7) * t391;
t672 = pkin(7) * t420;
t671 = pkin(7) * t429;
t670 = t587 * g(1);
t330 = t586 * t344 - t587 * t469;
t669 = qJ(2) * t330;
t348 = t586 * t377 - t587 * t431;
t668 = qJ(2) * t348;
t354 = t586 * t388 - t587 * t693;
t667 = qJ(2) * t354;
t585 = qJDD(1) * pkin(1);
t663 = t562 * t588;
t662 = t562 * t591;
t661 = t568 * t586;
t581 = t586 ^ 2;
t660 = t581 * t594;
t659 = t586 * t587;
t505 = -t586 * g(1) + t692 * t587;
t607 = -pkin(6) * t586 - t682;
t547 = t607 * qJD(1);
t478 = t547 * t633 + t505;
t560 = t593 * g(2) + t590 * g(3);
t536 = -t594 * qJ(2) + qJDD(2) + t560 - t585;
t509 = t607 * qJDD(1) + t536;
t496 = t592 * t509;
t620 = t569 * t636;
t641 = t592 * t594;
t413 = -t568 * pkin(3) - t535 * pkin(7) + t496 + (-pkin(3) * t581 * t641 + pkin(7) * t620 - t478) * t589;
t450 = t592 * t478 + t589 * t509;
t532 = -t569 * pkin(3) - pkin(7) * t618;
t567 = t687 * t660;
t414 = -pkin(3) * t567 - pkin(7) * t597 + t569 * t532 + t450;
t365 = -t591 * t413 + t588 * t414;
t621 = t468 * qJ(5) + t365;
t604 = -qJ(5) * t666 - t621;
t634 = qJD(5) * t528;
t333 = t604 - 0.2e1 * t634 + t691;
t658 = t588 * t333;
t608 = pkin(3) * t597 - pkin(7) * t567 + t670;
t611 = -t532 * t592 - t547;
t427 = (t543 + (t686 - t611) * qJD(1)) * t586 + t608;
t657 = t588 * t427;
t366 = t588 * t413 + t591 * t414;
t310 = -t591 * t365 + t588 * t366;
t654 = t589 * t310;
t477 = t670 + (t543 + (t686 + t547) * qJD(1)) * t586;
t653 = t589 * t477;
t622 = t589 * t641;
t556 = t581 * t622;
t533 = -t556 + t568;
t652 = t589 * t533;
t534 = -t556 - t568;
t651 = t589 * t534;
t650 = t590 * t536;
t649 = t591 * t333;
t648 = t591 * t427;
t645 = t592 * t310;
t644 = t592 * t477;
t643 = t592 * t533;
t642 = t592 * t534;
t640 = t593 * t536;
t639 = -pkin(1) * t342 + qJ(2) * t331;
t638 = -pkin(1) * t376 + qJ(2) * t349;
t637 = -pkin(1) * t386 + qJ(2) * t355;
t628 = t590 * qJDD(1);
t627 = t593 * qJDD(1);
t584 = t592 ^ 2;
t625 = t584 * t660;
t624 = t586 * t664;
t623 = t587 * t664;
t617 = -pkin(3) * t469 + pkin(7) * t393;
t616 = -pkin(3) * t431 + pkin(7) * t421;
t615 = -pkin(3) * t693 + pkin(7) * t430;
t605 = t593 * t594 + t628;
t614 = -pkin(5) * t605 + t593 * g(1);
t613 = -t536 + t585;
t311 = t588 * t365 + t591 * t366;
t449 = t589 * t478 - t496;
t504 = t692 * t586 + t670;
t461 = t586 * t504 + t587 * t505;
t516 = -t590 * t559 - t593 * t560;
t580 = t586 * t581;
t610 = t580 * t622;
t306 = -pkin(2) * t342 - pkin(3) * t391;
t606 = t587 * t556;
t398 = -t592 * t449 + t589 * t450;
t399 = t589 * t449 + t592 * t450;
t460 = t587 * t504 - t586 * t505;
t517 = t593 * t559 - t590 * t560;
t554 = -t590 * t594 + t627;
t582 = t587 ^ 2;
t546 = (t581 + t582) * t587 * t594;
t513 = -t590 * t546 + t587 * t627;
t515 = t593 * t546 + t587 * t628;
t603 = t467 * qJ(5) + 0.2e1 * qJD(5) * t526 + t562 * t494 + t366;
t601 = -pkin(1) * t330 + pkin(2) * t469 - pkin(6) * t344;
t600 = -pkin(1) * t348 + pkin(2) * t431 - pkin(6) * t377;
t599 = -pkin(1) * t354 + pkin(2) * t693 - pkin(6) * t388;
t372 = t427 + t688;
t575 = t582 * t594;
t574 = t582 * qJDD(1);
t573 = t581 * qJDD(1);
t552 = t575 - t660;
t551 = t575 + t660;
t550 = t574 - t573;
t549 = t574 + t573;
t545 = (t582 * t586 + t580) * t594;
t544 = t569 * t618;
t542 = t567 - t625;
t541 = t567 + t625;
t540 = t566 - t625;
t539 = -t567 - t566;
t538 = t567 - t566;
t537 = pkin(5) * t554 + t590 * g(1);
t523 = t554 * t659;
t522 = t605 * t659;
t520 = -t625 - t566;
t519 = 0.2e1 * t634;
t514 = -t593 * t545 - t586 * t628;
t512 = t590 * t545 - t586 * t627;
t511 = -t593 * t549 + t590 * t551;
t510 = t590 * t549 + t593 * t551;
t508 = (-t584 - t687) * t620;
t503 = -t525 + t557;
t502 = t524 - t557;
t499 = t544 - t597;
t498 = t544 + t597;
t497 = t609 - t535;
t492 = -t589 * t535 + t584 * t620;
t491 = (t687 * t569 * qJD(1) + t592 * t602) * t586;
t489 = (t630 + (qJD(3) - t569) * t635) * t589 * t586;
t487 = t592 * t539 - t651;
t486 = t592 * t538 + t652;
t485 = -t589 * t540 + t642;
t484 = t589 * t539 + t642;
t483 = -t589 * t538 + t643;
t482 = -t592 * t540 - t651;
t481 = -t525 + t524;
t480 = -t589 * t520 + t643;
t479 = t592 * t520 + t652;
t471 = t587 * t490 + t610;
t470 = t587 * t489 - t610;
t466 = (-t526 * t591 - t528 * t588) * t562;
t465 = (-t526 * t588 + t528 * t591) * t562;
t458 = t592 * t499 - t589 * t500;
t457 = -t589 * t497 - t592 * t498;
t456 = -t589 * t499 - t490;
t455 = t592 * t497 - t589 * t498;
t454 = t587 * t487 - t586 * t499;
t453 = t587 * t486 - t586 * t498;
t452 = t587 * t485 - t586 * t497;
t451 = t586 * t487 + t587 * t499;
t448 = t587 * t480 + t500 * t586;
t447 = t586 * t480 - t500 * t587;
t446 = t591 * t502 + t656;
t445 = -t588 * t503 + t646;
t444 = t588 * t502 - t647;
t443 = t591 * t503 + t655;
t442 = -t593 * t461 - t650;
t441 = t590 * t461 - t640;
t440 = -pkin(6) * t484 + t653;
t439 = t587 * t458 - t586 * t542;
t438 = t587 * t457 - t586 * t541;
t437 = t586 * t457 + t587 * t541;
t426 = -pkin(6) * t479 + t644;
t425 = t591 * t468 + t528 * t663;
t424 = t588 * t468 - t528 * t662;
t423 = -t588 * t467 + t526 * t662;
t422 = t591 * t467 + t526 * t663;
t418 = -pkin(2) * t484 + t449;
t417 = -t593 * t454 - t590 * t484;
t416 = t590 * t454 - t593 * t484;
t415 = -pkin(2) * t479 + t450;
t412 = -t589 * t465 + t592 * t466;
t411 = -t592 * t465 - t589 * t466;
t406 = -t593 * t448 - t590 * t479;
t405 = t590 * t448 - t593 * t479;
t404 = -pkin(4) * t693 + qJ(5) * t473;
t403 = t587 * t412 - t586 * t558;
t402 = t586 * t412 + t587 * t558;
t401 = -t593 * t438 - t590 * t455;
t400 = t590 * t438 - t593 * t455;
t397 = -t589 * t444 + t592 * t446;
t396 = -t589 * t443 + t592 * t445;
t395 = -t592 * t444 - t589 * t446;
t394 = -t592 * t443 - t589 * t445;
t392 = -t591 * t431 - t588 * t693;
t390 = -t588 * t431 + t591 * t693;
t387 = t648 - t671;
t384 = -pkin(1) * t451 - pkin(2) * t499 - pkin(6) * t487 + t644;
t383 = -t589 * t424 + t592 * t425;
t382 = -t589 * t422 + t592 * t423;
t381 = -t592 * t424 - t589 * t425;
t380 = -t592 * t422 - t589 * t423;
t379 = t657 - t672;
t378 = -pkin(1) * t447 + pkin(2) * t500 - pkin(6) * t480 - t653;
t375 = t587 * t399 + t586 * t477;
t374 = t586 * t399 - t587 * t477;
t371 = -pkin(6) * t455 - t398;
t370 = t587 * t383 - t624;
t369 = t587 * t382 + t624;
t368 = t586 * t383 + t623;
t367 = t586 * t382 - t623;
t363 = -qJ(5) * t488 + t372;
t362 = t587 * t397 - t586 * t432;
t361 = t587 * t396 - t586 * t436;
t360 = t586 * t397 + t587 * t432;
t359 = t586 * t396 + t587 * t436;
t358 = t615 + t657;
t357 = -t593 * t403 + t590 * t411;
t356 = t590 * t403 + t593 * t411;
t352 = -qJ(2) * t451 - t586 * t418 + t587 * t440;
t351 = t616 - t648;
t350 = -qJ(2) * t447 - t586 * t415 + t587 * t426;
t346 = -pkin(1) * t437 - pkin(2) * t541 - pkin(6) * t457 - t399;
t345 = -pkin(4) * t431 + qJ(5) * t472 - t608 + (t611 * qJD(1) - t692) * t586 - t688;
t343 = -t589 * t390 + t592 * t392;
t341 = -t592 * t390 - t589 * t392;
t339 = -t593 * t375 - t590 * t398;
t338 = t590 * t375 - t593 * t398;
t337 = -qJ(2) * t437 + t587 * t371 + t455 * t683;
t336 = -t524 * pkin(4) + t603;
t335 = t587 * t343 - t586 * t481;
t334 = t586 * t343 + t587 * t481;
t332 = -pkin(1) * t374 + pkin(2) * t477 - pkin(6) * t399;
t328 = -t593 * t370 + t590 * t381;
t327 = -t593 * t369 + t590 * t380;
t326 = t590 * t370 + t593 * t381;
t325 = t590 * t369 + t593 * t380;
t324 = -t593 * t362 + t590 * t395;
t323 = -t593 * t361 + t590 * t394;
t322 = t590 * t362 + t593 * t395;
t321 = t590 * t361 + t593 * t394;
t319 = t590 * t355 - t593 * t386;
t318 = pkin(5) * t319;
t317 = t591 * t363 - t588 * t404 - t671;
t316 = t519 + (-t436 + t666) * qJ(5) - t691 + t621;
t314 = t590 * t349 - t593 * t376;
t313 = -qJ(5) * t646 - t588 * t345 - t672;
t312 = pkin(5) * t314;
t309 = t366 + t690;
t308 = -qJ(2) * t374 + (-pkin(6) * t587 + t683) * t398;
t307 = -qJ(5) * t432 + (-t469 - t524) * pkin(4) + t603;
t305 = t365 + t689;
t304 = t588 * t363 + t591 * t404 + t615;
t303 = -qJ(5) * t655 + t591 * t345 + t616;
t302 = -pkin(3) * t427 + pkin(7) * t311;
t301 = -pkin(4) * t372 + qJ(5) * t336;
t300 = -pkin(4) * t436 + t306;
t299 = -t589 * t358 + t592 * t387 - t674;
t298 = -t589 * t351 + t592 * t379 - t675;
t297 = -t310 - t673;
t296 = (-t488 - t524) * pkin(4) + t603 + t690;
t295 = t311 + t617;
t294 = t519 - t604 + t689 - 0.2e1 * t691;
t293 = -t593 * t335 + t590 * t341;
t292 = t590 * t335 + t593 * t341;
t290 = t590 * t331 - t593 * t342;
t289 = pkin(5) * t290;
t288 = t591 * t336 - t658;
t287 = t588 * t336 + t649;
t286 = -t592 * t358 - t589 * t387 + t599;
t285 = t592 * t311 - t654;
t284 = t589 * t311 + t645;
t283 = -t592 * t351 - t589 * t379 + t600;
t282 = t587 * t285 + t586 * t427;
t281 = t586 * t285 - t587 * t427;
t280 = -t588 * t307 + t591 * t316 - t673;
t279 = -t589 * t304 + t592 * t317 - t674;
t278 = t591 * t307 + t588 * t316 + t617;
t277 = -t589 * t303 + t592 * t313 - t675;
t276 = t587 * t299 - t586 * t309 - t667;
t275 = -pkin(2) * t284 - pkin(3) * t310;
t274 = t587 * t298 - t586 * t305 - t668;
t273 = -t592 * t304 - t589 * t317 + t599;
t272 = -t589 * t287 + t592 * t288;
t271 = t592 * t287 + t589 * t288;
t270 = -t592 * t303 - t589 * t313 + t600;
t269 = -t589 * t295 + t592 * t297 - t676;
t268 = -pkin(7) * t287 - qJ(5) * t649 - t588 * t301;
t267 = t587 * t272 + t586 * t372;
t266 = t586 * t272 - t587 * t372;
t265 = -pkin(3) * t372 + pkin(7) * t288 - qJ(5) * t658 + t591 * t301;
t264 = -pkin(6) * t284 - pkin(7) * t645 - t589 * t302;
t263 = -t593 * t282 - t590 * t284;
t262 = t590 * t282 - t593 * t284;
t261 = t587 * t279 - t586 * t296 - t667;
t260 = -t592 * t295 - t589 * t297 + t601;
t259 = t587 * t277 - t586 * t294 - t668;
t258 = t587 * t269 - t586 * t306 - t669;
t257 = -pkin(2) * t271 - pkin(3) * t287 - pkin(4) * t333;
t256 = -t589 * t278 + t592 * t280 - t676;
t255 = -t592 * t278 - t589 * t280 + t601;
t254 = -pkin(1) * t281 + pkin(2) * t427 - pkin(6) * t285 + pkin(7) * t654 - t592 * t302;
t253 = -t593 * t267 - t590 * t271;
t252 = t590 * t267 - t593 * t271;
t251 = t587 * t256 - t586 * t300 - t669;
t250 = -qJ(2) * t281 + t587 * t264 - t586 * t275;
t249 = -pkin(6) * t271 - t589 * t265 + t592 * t268;
t248 = -pkin(1) * t266 + pkin(2) * t372 - pkin(6) * t272 - t592 * t265 - t589 * t268;
t247 = -qJ(2) * t266 + t587 * t249 - t586 * t257;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t460, 0, 0, 0, 0, 0, 0, t451, t447, t437, t374, 0, 0, 0, 0, 0, 0, t348, t354, t330, t281, 0, 0, 0, 0, 0, 0, t348, t354, t330, t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t554, -t605, 0, t516, 0, 0, 0, 0, 0, 0, t513, t512, t510, t441, 0, 0, 0, 0, 0, 0, t416, t405, t400, t338, 0, 0, 0, 0, 0, 0, t314, t319, t290, t262, 0, 0, 0, 0, 0, 0, t314, t319, t290, t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t605, t554, 0, t517, 0, 0, 0, 0, 0, 0, t515, t514, t511, t442, 0, 0, 0, 0, 0, 0, t417, t406, t401, t339, 0, 0, 0, 0, 0, 0, t315, t320, t291, t263, 0, 0, 0, 0, 0, 0, t315, t320, t291, t253; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t560, t559, 0, 0, t573, 0.2e1 * t586 * t629, 0, t574, 0, 0, -qJ(2) * t546 + t613 * t587, qJ(2) * t545 - t613 * t586, pkin(1) * t551 + qJ(2) * t549 + t461, -pkin(1) * t536 + qJ(2) * t461, t586 * t490 - t606, t586 * t458 + t587 * t542, t586 * t485 + t587 * t497, t586 * t489 + t606, t586 * t486 + t587 * t498, t587 * t568, -pkin(1) * t484 + qJ(2) * t454 + t587 * t418 + t586 * t440, -pkin(1) * t479 + qJ(2) * t448 + t587 * t415 + t586 * t426, qJ(2) * t438 + t586 * t371 + (-pkin(1) - t682) * t455, qJ(2) * t375 + (-pkin(1) + t607) * t398, t368, t334, t359, t367, t360, t402, t586 * t298 + t587 * t305 + t638, t586 * t299 + t587 * t309 + t637, t586 * t269 + t587 * t306 + t639, -pkin(1) * t284 + qJ(2) * t282 + t586 * t264 + t587 * t275, t368, t334, t359, t367, t360, t402, t586 * t277 + t587 * t294 + t638, t586 * t279 + t587 * t296 + t637, t586 * t256 + t587 * t300 + t639, -pkin(1) * t271 + qJ(2) * t267 + t586 * t249 + t587 * t257; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t605, 0, t554, 0, t614, -t537, -t517, -pkin(5) * t517, t522, t590 * t550 + t593 * t552, t512, -t522, -t513, 0, -pkin(5) * t515 + t593 * t504 + t586 * t650, -pkin(5) * t514 + t593 * t505 + t587 * t650, -pkin(5) * t511 + t590 * t460, -pkin(5) * t442 - (-pkin(1) * t593 - qJ(2) * t590) * t460, t590 * t471 + t593 * t492, t590 * t439 + t593 * t456, t590 * t452 + t593 * t482, t590 * t470 + t593 * t491, t590 * t453 + t593 * t483, t593 * t508 - t590 * t661, -pkin(5) * t417 + t590 * t352 + t593 * t384, -pkin(5) * t406 + t590 * t350 + t593 * t378, -pkin(5) * t401 + t590 * t337 + t593 * t346, -pkin(5) * t339 + t590 * t308 + t593 * t332, t326, t292, t321, t325, t322, t356, t590 * t274 + t593 * t283 - t678, t590 * t276 + t593 * t286 - t677, t590 * t258 + t593 * t260 - t679, -pkin(5) * t263 + t590 * t250 + t593 * t254, t326, t292, t321, t325, t322, t356, t590 * t259 + t593 * t270 - t678, t590 * t261 + t593 * t273 - t677, t590 * t251 + t593 * t255 - t679, -pkin(5) * t253 + t590 * t247 + t593 * t248; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t554, 0, t605, 0, t537, t614, t516, pkin(5) * t516, -t523, -t593 * t550 + t590 * t552, t514, t523, -t515, 0, pkin(5) * t513 + t590 * t504 - t586 * t640, pkin(5) * t512 + t590 * t505 - t587 * t640, pkin(5) * t510 - t593 * t460, pkin(5) * t441 - (-pkin(1) * t590 + qJ(2) * t593) * t460, -t593 * t471 + t590 * t492, -t593 * t439 + t590 * t456, -t593 * t452 + t590 * t482, -t593 * t470 + t590 * t491, -t593 * t453 + t590 * t483, t590 * t508 + t593 * t661, pkin(5) * t416 - t593 * t352 + t590 * t384, pkin(5) * t405 - t593 * t350 + t590 * t378, pkin(5) * t400 - t593 * t337 + t590 * t346, pkin(5) * t338 - t593 * t308 + t590 * t332, t328, t293, t323, t327, t324, t357, -t593 * t274 + t590 * t283 + t312, -t593 * t276 + t590 * t286 + t318, -t593 * t258 + t590 * t260 + t289, pkin(5) * t262 - t593 * t250 + t590 * t254, t328, t293, t323, t327, t324, t357, -t593 * t259 + t590 * t270 + t312, -t593 * t261 + t590 * t273 + t318, -t593 * t251 + t590 * t255 + t289, pkin(5) * t252 - t593 * t247 + t590 * t248;];
tauB_reg = t1;