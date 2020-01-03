% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPPR11_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:48:07
% EndTime: 2019-12-31 19:48:19
% DurationCPUTime: 9.87s
% Computational Cost: add. (28082->537), mult. (62512->751), div. (0->0), fcn. (37667->8), ass. (0->378)
t600 = sin(qJ(2));
t593 = t600 ^ 2;
t606 = qJD(1) ^ 2;
t585 = t593 * t606;
t605 = qJD(2) ^ 2;
t573 = -t585 - t605;
t603 = cos(qJ(2));
t635 = t603 * t606 * t600;
t568 = qJDD(2) - t635;
t658 = t603 * t568;
t512 = t600 * t573 + t658;
t645 = qJD(1) * qJD(2);
t582 = t603 * t645;
t643 = t600 * qJDD(1);
t556 = 0.2e1 * t582 + t643;
t601 = sin(qJ(1));
t604 = cos(qJ(1));
t458 = t601 * t512 + t604 * t556;
t714 = pkin(5) * t458;
t462 = t604 * t512 - t601 * t556;
t713 = pkin(5) * t462;
t594 = t603 ^ 2;
t587 = t594 * t606;
t575 = -t587 - t605;
t566 = qJDD(2) + t635;
t671 = t600 * t566;
t511 = -t603 * t575 + t671;
t631 = t600 * t645;
t641 = t603 * qJDD(1);
t559 = -0.2e1 * t631 + t641;
t457 = t601 * t511 - t604 * t559;
t712 = pkin(5) * t457;
t461 = t604 * t511 + t601 * t559;
t711 = pkin(5) * t461;
t659 = t603 * t566;
t502 = t600 * t575 + t659;
t710 = pkin(1) * t502;
t709 = pkin(6) * t502;
t708 = pkin(6) * t512;
t572 = -t585 + t605;
t506 = -t600 * t572 + t659;
t640 = t604 * qJDD(1);
t707 = t601 * t506 - t600 * t640;
t642 = t601 * qJDD(1);
t706 = t604 * t506 + t600 * t642;
t705 = 2 * qJD(3);
t704 = -2 * qJD(4);
t670 = t600 * t568;
t504 = -t603 * t573 + t670;
t703 = pkin(1) * t504;
t702 = pkin(6) * t504;
t701 = pkin(6) * t511;
t596 = sin(pkin(8));
t597 = cos(pkin(8));
t650 = qJD(1) * t603;
t542 = t596 * qJD(2) + t597 * t650;
t544 = t597 * qJD(2) - t596 * t650;
t495 = t544 * t542;
t557 = t582 + t643;
t692 = -t495 + t557;
t700 = t596 * t692;
t699 = t597 * t692;
t599 = sin(qJ(5));
t602 = cos(qJ(5));
t485 = t602 * t542 + t599 * t544;
t487 = -t599 * t542 + t602 * t544;
t429 = t487 * t485;
t545 = qJDD(5) + t557;
t693 = -t429 + t545;
t698 = t599 * t693;
t697 = t602 * t693;
t574 = t587 - t605;
t509 = -t603 * t574 + t670;
t696 = t601 * t509 + t603 * t640;
t695 = t604 * t509 - t601 * t641;
t617 = t557 + t582;
t694 = qJ(3) * t617;
t647 = t600 * qJD(1);
t691 = -pkin(2) * t631 + t647 * t705;
t558 = -t631 + t641;
t523 = t597 * qJDD(2) - t596 * t558;
t634 = t542 * t647;
t476 = t634 + t523;
t567 = pkin(3) * t647 - qJD(2) * qJ(4);
t569 = t601 * g(1) - t604 * g(2);
t614 = -qJDD(1) * pkin(1) - t569;
t608 = t614 - t691 - t694;
t685 = pkin(2) + qJ(4);
t402 = -t567 * t647 + (-pkin(3) * t594 - pkin(6)) * t606 - t685 * t558 + t608;
t570 = t604 * g(1) + t601 * g(2);
t536 = -t606 * pkin(1) + qJDD(1) * pkin(6) - t570;
t514 = t603 * g(3) + t600 * t536;
t684 = qJ(3) * t600;
t688 = pkin(2) * t603;
t619 = -t684 - t688;
t554 = t619 * qJD(1);
t613 = -qJDD(2) * pkin(2) - t605 * qJ(3) + t554 * t647 + qJDD(3) + t514;
t420 = -t566 * qJ(4) + (t557 - t582) * pkin(3) + t613;
t620 = -t596 * t402 + t597 * t420 + t544 * t704;
t356 = t597 * t402 + t596 * t420 + t542 * t704;
t690 = t600 * t574 + t658;
t589 = t600 * g(3);
t689 = -(qJD(1) * t554 + t536) * t603 + t605 * pkin(2) + t589;
t483 = t485 ^ 2;
t484 = t487 ^ 2;
t540 = t542 ^ 2;
t541 = t544 ^ 2;
t578 = qJD(5) + t647;
t576 = t578 ^ 2;
t651 = t593 + t594;
t561 = t651 * qJDD(1);
t564 = t585 + t587;
t493 = t601 * t561 + t604 * t564;
t687 = pkin(5) * t493;
t686 = t606 * pkin(6);
t683 = t578 * t599;
t682 = t578 * t602;
t331 = pkin(4) * t692 - t476 * pkin(7) + t620;
t522 = -t596 * qJDD(2) - t597 * t558;
t524 = pkin(4) * t647 - t544 * pkin(7);
t342 = -t540 * pkin(4) + t522 * pkin(7) - t524 * t647 + t356;
t299 = -t602 * t331 + t599 * t342;
t300 = t599 * t331 + t602 * t342;
t280 = -t602 * t299 + t599 * t300;
t681 = t596 * t280;
t609 = qJDD(2) * qJ(3) - t689;
t414 = qJDD(4) + t558 * pkin(3) - qJ(4) * t587 + (t705 + t567) * qJD(2) + t609;
t680 = t596 * t414;
t478 = t495 + t557;
t679 = t596 * t478;
t678 = t597 * t280;
t677 = t597 * t414;
t676 = t597 * t478;
t367 = -t522 * pkin(4) - t540 * pkin(7) + t544 * t524 + t414;
t675 = t599 * t367;
t416 = t429 + t545;
t674 = t599 * t416;
t535 = -t614 + t686;
t673 = t600 * t535;
t672 = t600 * t556;
t663 = t602 * t367;
t662 = t602 * t416;
t661 = t603 * t535;
t660 = t603 * t559;
t654 = -t541 - t585;
t653 = pkin(1) * t564 + pkin(6) * t561;
t646 = qJD(5) + t578;
t644 = qJD(3) * qJD(2);
t639 = t600 * t429;
t638 = t600 * t495;
t637 = t603 * t429;
t636 = t603 * t495;
t633 = t544 * t647;
t281 = t599 * t299 + t602 * t300;
t515 = t603 * t536 - t589;
t442 = t600 * t514 + t603 * t515;
t625 = -t602 * t522 + t599 * t523;
t498 = -t601 * t569 - t604 * t570;
t623 = t601 * t635;
t622 = t604 * t635;
t563 = -t601 * t606 + t640;
t621 = -pkin(5) * t563 - t601 * g(3);
t618 = pkin(2) * t600 - qJ(3) * t603;
t316 = t596 * t356 + t597 * t620;
t317 = t597 * t356 - t596 * t620;
t441 = t603 * t514 - t600 * t515;
t616 = t599 * t522 + t602 * t523;
t615 = t603 * t572 + t671;
t497 = t604 * t569 - t601 * t570;
t611 = t522 + t633;
t610 = (-qJD(5) + t578) * t487 - t625;
t411 = -t485 * qJD(5) + t616;
t443 = t609 + 0.2e1 * t644;
t607 = t558 * pkin(2) + t535 + t691;
t565 = -t585 + t587;
t562 = t604 * t606 + t642;
t552 = t618 * qJDD(1);
t547 = t651 * t645;
t538 = t603 * t557;
t537 = t600 * t557;
t529 = -pkin(5) * t562 + t604 * g(3);
t526 = -t541 + t585;
t525 = t540 - t585;
t521 = t601 * qJDD(2) + t604 * t547;
t520 = -t593 * t645 + t538;
t519 = -t604 * qJDD(2) + t601 * t547;
t518 = -t600 * t558 - t594 * t645;
t501 = t600 * t582 + t537;
t500 = (t558 - t631) * t603;
t494 = t604 * t561 - t601 * t564;
t492 = t541 - t540;
t491 = pkin(5) * t494;
t490 = t660 - t672;
t489 = t603 * t556 + t600 * t559;
t488 = -t585 - t540;
t482 = t604 * t520 - t623;
t481 = t604 * t518 + t623;
t480 = t601 * t520 + t622;
t479 = t601 * t518 - t622;
t475 = -t634 + t523;
t472 = -t522 + t633;
t471 = t578 * t485;
t470 = -t540 - t541;
t468 = (t542 * t597 - t544 * t596) * t647;
t467 = (t542 * t596 + t544 * t597) * t647;
t465 = -t484 + t576;
t464 = t483 - t576;
t455 = -t661 + t702;
t454 = -t673 - t709;
t453 = -t597 * t523 + t596 * t633;
t452 = -t596 * t523 - t597 * t633;
t451 = t596 * t522 - t597 * t634;
t450 = -t597 * t522 - t596 * t634;
t449 = t604 * t490 - t601 * t565;
t448 = t601 * t490 + t604 * t565;
t446 = t515 + t703;
t445 = t514 - t710;
t444 = -t484 - t576;
t439 = qJ(3) * t564 + t613;
t438 = -t600 * t467 + t538;
t437 = pkin(2) * t564 + t443;
t436 = t607 + t694;
t435 = -t597 * t525 + t679;
t434 = -t596 * t654 - t676;
t433 = t596 * t526 - t699;
t432 = -t596 * t525 - t676;
t431 = t597 * t654 - t679;
t430 = -t597 * t526 - t700;
t428 = t484 - t483;
t427 = -t686 + (-t558 - t559) * pkin(2) + t608;
t426 = (t556 + t617) * qJ(3) + t607;
t425 = t604 * t442 - t601 * t535;
t424 = t601 * t442 + t604 * t535;
t423 = -t576 - t483;
t422 = t597 * t488 - t700;
t421 = t596 * t488 + t699;
t419 = -t600 * t452 + t636;
t418 = -t600 * t450 - t636;
t410 = -t487 * qJD(5) - t625;
t409 = t596 * t476 + t597 * t611;
t408 = t597 * t472 + t596 * t475;
t407 = -t597 * t476 + t596 * t611;
t406 = t596 * t472 - t597 * t475;
t405 = (-t485 * t602 + t487 * t599) * t578;
t404 = (-t485 * t599 - t487 * t602) * t578;
t403 = pkin(2) * t566 + qJ(3) * t575 - t613 + t710;
t401 = -t703 + pkin(2) * t573 - 0.2e1 * t644 + (-qJDD(2) - t568) * qJ(3) + t689;
t399 = -t483 - t484;
t398 = t603 * t443 + t600 * t613;
t397 = t600 * t443 - t603 * t613;
t396 = -t600 * t430 + t603 * t476;
t395 = t600 * t431 + t603 * t475;
t394 = -t600 * t432 + t603 * t611;
t393 = -t603 * t431 + t600 * t475;
t392 = -pkin(2) * t672 + t603 * t426 - t702;
t391 = -qJ(3) * t660 - t600 * t427 + t709;
t390 = -t600 * t406 + t603 * t492;
t389 = t600 * t421 + t603 * t472;
t388 = -t603 * t421 + t600 * t472;
t387 = t411 + t471;
t386 = t411 - t471;
t385 = -t646 * t485 + t616;
t382 = t646 * t487 + t625;
t381 = -t600 * t437 + t603 * t439;
t380 = t602 * t464 - t674;
t379 = -t599 * t465 + t697;
t378 = t599 * t464 + t662;
t377 = t602 * t465 + t698;
t376 = t602 * t411 - t487 * t683;
t375 = t599 * t411 + t487 * t682;
t374 = -t599 * t410 + t485 * t682;
t373 = t602 * t410 + t485 * t683;
t372 = t600 * t407 + t603 * t470;
t371 = -t603 * t407 + t600 * t470;
t370 = -t599 * t444 - t662;
t369 = t602 * t444 - t674;
t366 = t602 * t423 - t698;
t365 = t599 * t423 + t697;
t363 = t604 * t398 - t601 * t436;
t362 = t601 * t398 + t604 * t436;
t361 = pkin(3) * t407 - qJ(3) * t409;
t360 = t596 * t404 - t597 * t405;
t359 = -t597 * t404 - t596 * t405;
t358 = t604 * t395 + t601 * t434;
t357 = t601 * t395 - t604 * t434;
t353 = -t600 * t359 + t603 * t545;
t352 = -pkin(1) * t397 + pkin(2) * t613 - qJ(3) * t443;
t351 = t604 * t389 + t601 * t422;
t350 = t601 * t389 - t604 * t422;
t349 = t604 * t372 + t601 * t409;
t348 = t601 * t372 - t604 * t409;
t347 = -pkin(6) * t397 - t436 * t618;
t346 = t599 * t387 + t602 * t610;
t345 = -t602 * t382 - t599 * t386;
t344 = -t602 * t387 + t599 * t610;
t343 = -t599 * t382 + t602 * t386;
t341 = pkin(3) * t475 - t685 * t434 - t680;
t339 = t596 * t378 - t597 * t380;
t338 = t596 * t377 - t597 * t379;
t337 = -t597 * t378 - t596 * t380;
t336 = -t597 * t377 - t596 * t379;
t335 = t596 * t375 - t597 * t376;
t334 = t596 * t373 - t597 * t374;
t333 = -t597 * t375 - t596 * t376;
t332 = -t597 * t373 - t596 * t374;
t329 = -t596 * t369 + t597 * t370;
t328 = t597 * t369 + t596 * t370;
t327 = pkin(3) * t472 - t685 * t422 + t677;
t326 = -pkin(7) * t369 + t663;
t325 = -pkin(7) * t365 + t675;
t324 = pkin(3) * t431 - qJ(3) * t434 - t356;
t323 = -t596 * t365 + t597 * t366;
t322 = t597 * t365 + t596 * t366;
t321 = -t600 * t333 + t637;
t320 = -t600 * t332 - t637;
t319 = pkin(3) * t421 - qJ(3) * t422 + t620;
t318 = -pkin(1) * t393 - qJ(3) * t475 + t685 * t431 - t677;
t315 = -t600 * t336 + t603 * t387;
t314 = -t600 * t337 + t603 * t610;
t313 = -pkin(1) * t388 - qJ(3) * t472 + t685 * t421 - t680;
t312 = t600 * t328 + t603 * t385;
t311 = -t603 * t328 + t600 * t385;
t310 = -pkin(4) * t385 + pkin(7) * t370 + t675;
t309 = -pkin(4) * t382 + pkin(7) * t366 - t663;
t308 = t600 * t322 + t603 * t382;
t307 = -t603 * t322 + t600 * t382;
t306 = t600 * t316 + t603 * t414;
t305 = -t603 * t316 + t600 * t414;
t304 = -t596 * t344 + t597 * t346;
t303 = t596 * t343 - t597 * t345;
t302 = t597 * t344 + t596 * t346;
t301 = -t597 * t343 - t596 * t345;
t297 = -t600 * t301 + t603 * t428;
t296 = pkin(3) * t470 - t685 * t409 - t317;
t295 = t600 * t302 + t603 * t399;
t294 = -t603 * t302 + t600 * t399;
t293 = t604 * t312 + t601 * t329;
t292 = t601 * t312 - t604 * t329;
t291 = -pkin(6) * t393 + t603 * t324 - t600 * t341;
t290 = -pkin(1) * t371 - qJ(3) * t470 + t685 * t407 + t316;
t289 = -pkin(6) * t388 + t603 * t319 - t600 * t327;
t288 = pkin(3) * t316 - qJ(3) * t317;
t287 = t604 * t308 + t601 * t323;
t286 = t601 * t308 - t604 * t323;
t285 = pkin(3) * t414 - t685 * t317;
t284 = t604 * t306 + t601 * t317;
t283 = t601 * t306 - t604 * t317;
t282 = -pkin(6) * t371 - t600 * t296 + t603 * t361;
t278 = t604 * t295 + t601 * t304;
t277 = t601 * t295 - t604 * t304;
t276 = pkin(3) * t328 + pkin(4) * t369 - qJ(3) * t329 - t300;
t275 = -pkin(4) * t367 + pkin(7) * t281;
t274 = pkin(3) * t302 + pkin(4) * t344 - qJ(3) * t304;
t273 = -pkin(7) * t344 - t280;
t272 = -pkin(1) * t305 - qJ(3) * t414 + t685 * t316;
t271 = pkin(3) * t322 + pkin(4) * t365 - qJ(3) * t323 - t299;
t270 = pkin(3) * t385 - t597 * t310 - t596 * t326 - t685 * t329;
t269 = -pkin(4) * t399 + pkin(7) * t346 + t281;
t268 = pkin(3) * t382 - t597 * t309 - t685 * t323 - t596 * t325;
t267 = -pkin(1) * t311 - qJ(3) * t385 + t596 * t310 - t597 * t326 + t685 * t328;
t266 = -pkin(1) * t307 - qJ(3) * t382 + t596 * t309 + t685 * t322 - t597 * t325;
t265 = -pkin(6) * t305 - t600 * t285 + t603 * t288;
t264 = t597 * t281 - t681;
t263 = t596 * t281 + t678;
t262 = t600 * t263 + t603 * t367;
t261 = -t603 * t263 + t600 * t367;
t260 = -pkin(6) * t311 - t600 * t270 + t603 * t276;
t259 = -pkin(6) * t307 - t600 * t268 + t603 * t271;
t258 = pkin(3) * t399 - t597 * t269 - t596 * t273 - t685 * t304;
t257 = t604 * t262 + t601 * t264;
t256 = t601 * t262 - t604 * t264;
t255 = -pkin(1) * t294 - qJ(3) * t399 + t596 * t269 - t597 * t273 + t685 * t302;
t254 = pkin(3) * t263 + pkin(4) * t280 - qJ(3) * t264;
t253 = -pkin(6) * t294 - t600 * t258 + t603 * t274;
t252 = pkin(3) * t367 + pkin(7) * t681 - t685 * t264 - t597 * t275;
t251 = -pkin(1) * t261 + pkin(7) * t678 - qJ(3) * t367 + t685 * t263 + t596 * t275;
t250 = -pkin(6) * t261 - t600 * t252 + t603 * t254;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t562, -t563, 0, t498, 0, 0, 0, 0, 0, 0, -t461, -t462, t494, t425, 0, 0, 0, 0, 0, 0, t494, t461, t462, t363, 0, 0, 0, 0, 0, 0, t351, t358, t349, t284, 0, 0, 0, 0, 0, 0, t287, t293, t278, t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t563, -t562, 0, t497, 0, 0, 0, 0, 0, 0, -t457, -t458, t493, t424, 0, 0, 0, 0, 0, 0, t493, t457, t458, t362, 0, 0, 0, 0, 0, 0, t350, t357, t348, t283, 0, 0, 0, 0, 0, 0, t286, t292, t277, t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t502, -t504, 0, -t441, 0, 0, 0, 0, 0, 0, 0, -t502, t504, t397, 0, 0, 0, 0, 0, 0, t388, t393, t371, t305, 0, 0, 0, 0, 0, 0, t307, t311, t294, t261; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t563, 0, -t562, 0, t621, -t529, -t497, -pkin(5) * t497, t482, t449, t706, t481, -t695, t521, -t601 * t445 + t604 * t454 + t712, -t601 * t446 + t604 * t455 + t714, t604 * t441 - t687, -pkin(5) * t424 - (pkin(1) * t601 - pkin(6) * t604) * t441, t521, -t706, t695, t482, t449, t481, t604 * t381 - t601 * t552 - t687, t604 * t391 - t601 * t403 - t712, t604 * t392 - t601 * t401 - t714, -pkin(5) * t362 + t604 * t347 - t601 * t352, t604 * t419 - t601 * t453, t604 * t390 - t601 * t408, t604 * t396 - t601 * t433, t604 * t418 - t601 * t451, t604 * t394 - t601 * t435, t604 * t438 - t601 * t468, -pkin(5) * t350 + t604 * t289 - t601 * t313, -pkin(5) * t357 + t604 * t291 - t601 * t318, -pkin(5) * t348 + t604 * t282 - t601 * t290, -pkin(5) * t283 + t604 * t265 - t601 * t272, t604 * t321 - t601 * t335, t604 * t297 - t601 * t303, t604 * t315 - t601 * t338, t604 * t320 - t601 * t334, t604 * t314 - t601 * t339, t604 * t353 - t601 * t360, -pkin(5) * t286 + t604 * t259 - t601 * t266, -pkin(5) * t292 + t604 * t260 - t601 * t267, -pkin(5) * t277 + t604 * t253 - t601 * t255, -pkin(5) * t256 + t604 * t250 - t601 * t251; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t562, 0, t563, 0, t529, t621, t498, pkin(5) * t498, t480, t448, t707, t479, -t696, t519, t604 * t445 + t601 * t454 - t711, t604 * t446 + t601 * t455 - t713, t601 * t441 + t491, pkin(5) * t425 - (-pkin(1) * t604 - pkin(6) * t601) * t441, t519, -t707, t696, t480, t448, t479, t601 * t381 + t604 * t552 + t491, t601 * t391 + t604 * t403 + t711, t601 * t392 + t604 * t401 + t713, pkin(5) * t363 + t601 * t347 + t604 * t352, t601 * t419 + t604 * t453, t601 * t390 + t604 * t408, t601 * t396 + t604 * t433, t601 * t418 + t604 * t451, t601 * t394 + t604 * t435, t601 * t438 + t604 * t468, pkin(5) * t351 + t601 * t289 + t604 * t313, pkin(5) * t358 + t601 * t291 + t604 * t318, pkin(5) * t349 + t601 * t282 + t604 * t290, pkin(5) * t284 + t601 * t265 + t604 * t272, t601 * t321 + t604 * t335, t601 * t297 + t604 * t303, t601 * t315 + t604 * t338, t601 * t320 + t604 * t334, t601 * t314 + t604 * t339, t601 * t353 + t604 * t360, pkin(5) * t287 + t601 * t259 + t604 * t266, pkin(5) * t293 + t601 * t260 + t604 * t267, pkin(5) * t278 + t601 * t253 + t604 * t255, pkin(5) * t257 + t601 * t250 + t604 * t251; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t569, t570, 0, 0, t501, t489, t615, t500, t690, 0, pkin(1) * t559 + t661 - t701, -pkin(1) * t556 - t673 - t708, t442 + t653, pkin(1) * t535 + pkin(6) * t442, 0, -t615, -t690, t501, t489, t500, t603 * t437 + t600 * t439 + t653, t701 + t603 * t427 + (-pkin(1) - t684) * t559, t708 + t600 * t426 + (pkin(1) + t688) * t556, pkin(6) * t398 + (pkin(1) - t619) * t436, t603 * t452 + t638, t603 * t406 + t600 * t492, t603 * t430 + t600 * t476, t603 * t450 - t638, t603 * t432 + t600 * t611, t603 * t467 + t537, -pkin(1) * t422 + pkin(6) * t389 + t600 * t319 + t603 * t327, -pkin(1) * t434 + pkin(6) * t395 + t600 * t324 + t603 * t341, -pkin(1) * t409 + pkin(6) * t372 + t603 * t296 + t600 * t361, -pkin(1) * t317 + pkin(6) * t306 + t603 * t285 + t600 * t288, t603 * t333 + t639, t603 * t301 + t600 * t428, t603 * t336 + t600 * t387, t603 * t332 - t639, t603 * t337 + t600 * t610, t603 * t359 + t600 * t545, -pkin(1) * t323 + pkin(6) * t308 + t603 * t268 + t600 * t271, -pkin(1) * t329 + pkin(6) * t312 + t603 * t270 + t600 * t276, -pkin(1) * t304 + pkin(6) * t295 + t603 * t258 + t600 * t274, -pkin(1) * t264 + pkin(6) * t262 + t603 * t252 + t600 * t254;];
tauB_reg = t1;