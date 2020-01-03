% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRP11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRP11_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:11
% EndTime: 2019-12-31 20:14:22
% DurationCPUTime: 8.56s
% Computational Cost: add. (12623->475), mult. (26366->600), div. (0->0), fcn. (15062->6), ass. (0->358)
t597 = sin(qJ(2));
t653 = t597 * qJD(1);
t578 = qJD(4) + t653;
t575 = t578 ^ 2;
t596 = sin(qJ(4));
t599 = cos(qJ(4));
t600 = cos(qJ(2));
t654 = qJD(1) * t600;
t547 = t596 * qJD(2) + t599 * t654;
t697 = t547 ^ 2;
t519 = t697 - t575;
t549 = t599 * qJD(2) - t596 * t654;
t494 = t549 * t547;
t650 = qJD(1) * qJD(2);
t582 = t600 * t650;
t648 = t597 * qJDD(1);
t555 = t582 + t648;
t540 = qJDD(4) + t555;
t708 = t494 + t540;
t668 = t599 * t708;
t411 = t596 * t519 + t668;
t524 = t578 * t549;
t637 = t597 * t650;
t646 = t600 * qJDD(1);
t556 = -t637 + t646;
t630 = t596 * qJDD(2) + t599 * t556;
t613 = t549 * qJD(4) + t630;
t431 = -t524 + t613;
t367 = t597 * t411 - t600 * t431;
t684 = t596 * t708;
t415 = t599 * t519 - t684;
t598 = sin(qJ(1));
t601 = cos(qJ(1));
t769 = t598 * t367 - t601 * t415;
t768 = t601 * t367 + t598 * t415;
t539 = t549 ^ 2;
t643 = t575 + t539;
t406 = t596 * t643 - t668;
t767 = pkin(1) * t406;
t592 = t597 ^ 2;
t603 = qJD(1) ^ 2;
t584 = t592 * t603;
t602 = qJD(2) ^ 2;
t572 = -t584 - t602;
t639 = t600 * t603 * t597;
t565 = qJDD(2) - t639;
t661 = t600 * t565;
t510 = t597 * t572 + t661;
t554 = 0.2e1 * t582 + t648;
t456 = t598 * t510 + t601 * t554;
t766 = pkin(5) * t456;
t460 = t601 * t510 - t598 * t554;
t765 = pkin(5) * t460;
t394 = t599 * t643 + t684;
t764 = t597 * t394;
t763 = t598 * t406;
t761 = t600 * t394;
t760 = t601 * t406;
t695 = pkin(2) + pkin(7);
t758 = t695 * t394;
t757 = t695 * t406;
t756 = t600 * t411 + t597 * t431;
t755 = -pkin(3) * t394 - qJ(3) * t406;
t593 = t600 ^ 2;
t586 = t593 * t603;
t574 = -t586 - t602;
t564 = qJDD(2) + t639;
t680 = t597 * t564;
t509 = -t600 * t574 + t680;
t557 = -0.2e1 * t637 + t646;
t455 = t598 * t509 - t601 * t557;
t754 = pkin(5) * t455;
t459 = t601 * t509 + t598 * t557;
t753 = pkin(5) * t459;
t709 = -t494 + t540;
t463 = t596 * t709;
t521 = -t539 + t575;
t731 = -t599 * t521 - t463;
t752 = t597 * t731;
t667 = t599 * t709;
t730 = t596 * t521 - t667;
t751 = t598 * t730;
t750 = t600 * t731;
t749 = t601 * t730;
t662 = t600 * t564;
t500 = t597 * t574 + t662;
t748 = pkin(1) * t500;
t705 = -t697 - t575;
t722 = t599 * t705 - t463;
t747 = pkin(1) * t722;
t464 = -t697 - t539;
t746 = pkin(3) * t464;
t745 = pkin(6) * t500;
t744 = pkin(6) * t510;
t743 = qJ(3) * t464;
t618 = t599 * qJDD(2) - t596 * t556;
t480 = -t547 * qJD(4) + t618;
t523 = t547 * t578;
t435 = t480 - t523;
t742 = qJ(5) * t435;
t741 = t597 * t464;
t490 = t539 - t697;
t740 = t597 * t490;
t719 = t596 * t705 + t667;
t739 = t597 * t719;
t738 = t598 * t722;
t737 = t600 * t464;
t736 = t600 * t490;
t735 = t600 * t719;
t734 = t601 * t722;
t733 = t695 * t719;
t732 = t695 * t722;
t571 = -t584 + t602;
t504 = -t597 * t571 + t662;
t645 = t601 * qJDD(1);
t729 = t598 * t504 - t597 * t645;
t647 = t598 * qJDD(1);
t728 = t601 * t504 + t597 * t647;
t727 = pkin(3) * t719 - qJ(3) * t722;
t726 = 2 * qJD(3);
t679 = t597 * t565;
t502 = -t600 * t572 + t679;
t725 = pkin(1) * t502;
t724 = pkin(6) * t502;
t723 = pkin(6) * t509;
t573 = t586 - t602;
t507 = -t600 * t573 + t679;
t721 = t598 * t507 + t600 * t645;
t720 = t601 * t507 - t598 * t646;
t687 = t578 * t599;
t625 = -t547 * t687 - t596 * t613;
t688 = t578 * t596;
t642 = t547 * t688;
t615 = t599 * t613 - t642;
t640 = t600 * t494;
t700 = -t597 * t615 - t640;
t718 = t598 * t700 + t601 * t625;
t612 = (t547 * t599 - t549 * t596) * t578;
t514 = t549 * t687;
t624 = t514 + t642;
t701 = t600 * t540 - t597 * t624;
t717 = t598 * t701 + t601 * t612;
t716 = -t598 * t625 + t601 * t700;
t715 = -t598 * t612 + t601 * t701;
t621 = t555 + t582;
t710 = t621 * qJ(3);
t707 = -pkin(2) * t637 + t653 * t726;
t706 = t597 * t573 + t661;
t704 = pkin(4) * t613 - t742;
t426 = -t599 * t480 + t549 * t688;
t425 = -t596 * t480 - t514;
t616 = -t597 * t425 + t640;
t703 = t601 * t426 + t598 * t616;
t702 = t597 * t540 + t600 * t624;
t641 = t597 * t494;
t699 = t600 * t615 - t641;
t698 = -t598 * t426 + t601 * t616;
t696 = 2 * qJD(5);
t694 = pkin(2) * t600;
t693 = pkin(4) * t596;
t692 = pkin(4) * t599;
t655 = t592 + t593;
t559 = t655 * qJDD(1);
t562 = t584 + t586;
t488 = t598 * t559 + t601 * t562;
t691 = pkin(5) * t488;
t690 = t603 * pkin(6);
t689 = qJ(3) * t597;
t567 = pkin(3) * t653 - qJD(2) * pkin(7);
t623 = -t689 - t694;
t552 = t623 * qJD(1);
t569 = t601 * g(1) + t598 * g(2);
t533 = -t603 * pkin(1) + qJDD(1) * pkin(6) - t569;
t656 = t597 * g(3) - t600 * t533;
t611 = t602 * pkin(2) - t552 * t654 + t656;
t609 = qJDD(2) * qJ(3) - t611;
t607 = t556 * pkin(3) - pkin(7) * t586 + t609;
t398 = (t726 + t567) * qJD(2) + t607;
t686 = t596 * t398;
t430 = t524 + t613;
t685 = t596 * t430;
t568 = t598 * g(1) - t601 * g(2);
t617 = -qJDD(1) * pkin(1) - t568;
t532 = -t617 + t690;
t682 = t597 * t532;
t681 = t597 * t554;
t671 = t599 * t398;
t670 = t599 * t430;
t652 = -qJD(4) + t578;
t438 = t652 * t547 + t618;
t669 = t599 * t438;
t664 = t600 * t532;
t663 = t600 * t557;
t606 = t617 - t707 - t710;
t386 = -t567 * t653 + (-pkin(3) * t593 - pkin(6)) * t603 - t695 * t556 + t606;
t512 = t600 * g(3) + t597 * t533;
t614 = -qJDD(2) * pkin(2) - t602 * qJ(3) + t552 * t653 + qJDD(3) + t512;
t399 = -t564 * pkin(7) + (t555 - t582) * pkin(3) + t614;
t349 = t599 * t386 + t596 * t399;
t658 = -t596 * t386 + t599 * t399;
t657 = pkin(1) * t562 + pkin(6) * t559;
t651 = qJD(4) + t578;
t649 = qJD(3) * qJD(2);
t632 = qJ(5) * t596 + pkin(3);
t631 = qJ(5) * t599 - qJ(3);
t441 = t597 * t512 - t600 * t656;
t496 = -t598 * t568 - t601 * t569;
t629 = t598 * t639;
t628 = t601 * t639;
t561 = -t598 * t603 + t645;
t627 = -pkin(5) * t561 - t598 * g(3);
t626 = t600 * t425 + t641;
t622 = pkin(2) * t597 - qJ(3) * t600;
t485 = t547 * pkin(4) - t549 * qJ(5);
t620 = t540 * qJ(5) - t547 * t485 + t578 * t696 + t349;
t318 = t596 * t349 + t599 * t658;
t319 = t599 * t349 - t596 * t658;
t440 = t600 * t512 + t597 * t656;
t619 = t600 * t571 + t680;
t495 = t601 * t568 - t598 * t569;
t610 = -t540 * pkin(4) - t575 * qJ(5) + t549 * t485 + qJDD(5) - t658;
t608 = -t651 * t547 + t618;
t442 = t609 + 0.2e1 * t649;
t605 = t556 * pkin(2) + t532 + t707;
t587 = -0.2e1 * t649;
t604 = -qJD(2) * t567 + t549 * t696 + t587 - t607 - t704;
t563 = t586 - t584;
t560 = t601 * t603 + t647;
t550 = t622 * qJDD(1);
t542 = t655 * t650;
t528 = -pkin(5) * t560 + t601 * g(3);
t518 = t598 * qJDD(2) + t601 * t542;
t517 = t600 * t555 - t592 * t650;
t516 = -t601 * qJDD(2) + t598 * t542;
t515 = -t597 * t556 - t593 * t650;
t499 = t621 * t597;
t498 = (t556 - t637) * t600;
t489 = t601 * t559 - t598 * t562;
t486 = pkin(5) * t489;
t484 = t663 - t681;
t483 = t600 * t554 + t597 * t557;
t478 = t601 * t517 - t629;
t477 = t601 * t515 + t629;
t476 = t598 * t517 + t628;
t475 = t598 * t515 - t628;
t453 = -t664 + t724;
t452 = -t682 - t745;
t447 = t601 * t484 - t598 * t563;
t446 = t598 * t484 + t601 * t563;
t444 = -t656 + t725;
t443 = t512 - t748;
t434 = t480 + t523;
t433 = t652 * t549 - t630;
t432 = t651 * t549 + t630;
t429 = t596 * t438;
t428 = qJ(3) * t562 + t614;
t427 = pkin(2) * t562 + t442;
t420 = t605 + t710;
t403 = -t690 + (-t556 - t557) * pkin(2) + t606;
t402 = (t554 + t621) * qJ(3) + t605;
t401 = t601 * t441 - t598 * t532;
t400 = t598 * t441 + t601 * t532;
t391 = pkin(2) * t564 + qJ(3) * t574 - t614 + t748;
t385 = -t725 + pkin(2) * t572 + t587 + (-qJDD(2) - t565) * qJ(3) + t611;
t382 = t600 * t442 + t597 * t614;
t381 = t597 * t442 - t600 * t614;
t380 = t599 * t433 + t429;
t379 = t596 * t608 + t670;
t378 = -t599 * t431 + t429;
t377 = -t596 * t435 - t670;
t376 = t596 * t433 - t669;
t375 = -t599 * t608 + t685;
t374 = -t596 * t431 - t669;
t373 = t599 * t435 - t685;
t372 = -pkin(2) * t681 + t600 * t402 - t724;
t371 = -qJ(3) * t663 - t597 * t403 + t745;
t370 = -t597 * t427 + t600 * t428;
t369 = t600 * t438 - t752;
t368 = t600 * t434 - t752;
t365 = t600 * t608 - t764;
t364 = t600 * t432 + t739;
t363 = t597 * t608 + t761;
t362 = t597 * t432 - t735;
t361 = -t600 * t435 + t764;
t360 = t600 * t430 + t739;
t359 = -t597 * t435 - t761;
t358 = t597 * t430 - t735;
t357 = -t597 * t373 - t736;
t356 = -t597 * t375 + t736;
t355 = t597 * t374 + t737;
t354 = t597 * t376 + t737;
t353 = -t600 * t374 + t741;
t352 = -t600 * t376 + t741;
t351 = t601 * t382 - t598 * t420;
t350 = t598 * t382 + t601 * t420;
t346 = -pkin(1) * t381 + pkin(2) * t614 - qJ(3) * t442;
t345 = (pkin(4) * t578 - (2 * qJD(5))) * t549 + t398 + t704;
t344 = t601 * t365 + t763;
t343 = t601 * t364 + t738;
t342 = t598 * t365 - t760;
t341 = t598 * t364 - t734;
t340 = t601 * t361 - t763;
t339 = t601 * t360 + t738;
t338 = t598 * t361 + t760;
t337 = t598 * t360 - t734;
t336 = pkin(3) * t376 - qJ(3) * t380;
t335 = -pkin(6) * t381 - t622 * t420;
t333 = -t575 * pkin(4) + t620;
t332 = (-t432 - t524) * pkin(4) + t604;
t331 = -pkin(4) * t524 + t604 + t742;
t330 = t601 * t355 + t598 * t378;
t329 = t601 * t354 + t598 * t380;
t328 = t598 * t355 - t601 * t378;
t327 = t598 * t354 - t601 * t380;
t326 = -qJ(5) * t464 + t610;
t325 = pkin(3) * t608 - t686 - t757;
t324 = (-t464 - t575) * pkin(4) + t620;
t323 = pkin(3) * t430 + t671 - t732;
t322 = -t349 + t755;
t321 = t658 + t727;
t320 = pkin(3) * t374 - pkin(4) * t438 - qJ(3) * t378 - qJ(5) * t431;
t317 = -pkin(1) * t363 - qJ(3) * t608 - t671 - t758;
t316 = -pkin(1) * t358 - qJ(3) * t430 - t686 + t733;
t315 = t597 * t318 + t600 * t398;
t314 = -t600 * t318 + t597 * t398;
t313 = pkin(4) * t709 + qJ(5) * t705 - t610 + t727;
t312 = t539 * pkin(4) + qJ(5) * t708 + t620 - t755;
t311 = t599 * t333 + t596 * t610;
t310 = t596 * t333 - t599 * t610;
t309 = -t599 * t332 + t632 * t432 - t732;
t308 = -t596 * t331 + (-pkin(3) - t692) * t435 + t757;
t307 = -t695 * t380 - t319 + t746;
t306 = -pkin(1) * t362 + t596 * t332 + t631 * t432 + t733;
t305 = -pkin(1) * t359 - t599 * t331 + (qJ(3) + t693) * t435 + t758;
t304 = t597 * t310 + t600 * t345;
t303 = -t600 * t310 + t597 * t345;
t302 = pkin(3) * t318 - qJ(3) * t319;
t301 = -pkin(6) * t363 + t600 * t322 - t597 * t325;
t300 = -pkin(6) * t358 + t600 * t321 - t597 * t323;
t299 = -pkin(1) * t352 + t695 * t376 + t318 - t743;
t298 = pkin(3) * t398 - t695 * t319;
t297 = t601 * t315 + t598 * t319;
t296 = t598 * t315 - t601 * t319;
t295 = -t599 * t324 - t596 * t326 - t695 * t378 + t746;
t294 = -pkin(1) * t353 + t596 * t324 - t599 * t326 + t695 * t374 - t743;
t293 = -pkin(6) * t352 - t597 * t307 + t600 * t336;
t292 = -pkin(6) * t362 - t597 * t309 + t600 * t313;
t291 = -pkin(6) * t359 - t597 * t308 + t600 * t312;
t290 = t601 * t304 + t598 * t311;
t289 = t598 * t304 - t601 * t311;
t288 = -pkin(1) * t314 - qJ(3) * t398 + t695 * t318;
t287 = -pkin(6) * t353 - t597 * t295 + t600 * t320;
t286 = pkin(3) * t310 - pkin(4) * t610 - qJ(3) * t311 + qJ(5) * t333;
t285 = -t695 * t311 + (t632 + t692) * t345;
t284 = -pkin(6) * t314 - t597 * t298 + t600 * t302;
t283 = -pkin(1) * t303 + t695 * t310 + (t631 - t693) * t345;
t282 = -pkin(6) * t303 - t597 * t285 + t600 * t286;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t560, -t561, 0, t496, 0, 0, 0, 0, 0, 0, -t459, -t460, t489, t401, 0, 0, 0, 0, 0, 0, t489, t459, t460, t351, 0, 0, 0, 0, 0, 0, t339, t344, t329, t297, 0, 0, 0, 0, 0, 0, t343, t330, t340, t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t561, -t560, 0, t495, 0, 0, 0, 0, 0, 0, -t455, -t456, t488, t400, 0, 0, 0, 0, 0, 0, t488, t455, t456, t350, 0, 0, 0, 0, 0, 0, t337, t342, t327, t296, 0, 0, 0, 0, 0, 0, t341, t328, t338, t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t500, -t502, 0, -t440, 0, 0, 0, 0, 0, 0, 0, -t500, t502, t381, 0, 0, 0, 0, 0, 0, t358, t363, t352, t314, 0, 0, 0, 0, 0, 0, t362, t353, t359, t303; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t561, 0, -t560, 0, t627, -t528, -t495, -pkin(5) * t495, t478, t447, t728, t477, -t720, t518, -t598 * t443 + t601 * t452 + t754, -t598 * t444 + t601 * t453 + t766, t601 * t440 - t691, -pkin(5) * t400 - (pkin(1) * t598 - pkin(6) * t601) * t440, t518, -t728, t720, t478, t447, t477, t601 * t370 - t598 * t550 - t691, t601 * t371 - t598 * t391 - t754, t601 * t372 - t598 * t385 - t766, -pkin(5) * t350 + t601 * t335 - t598 * t346, t698, t601 * t356 - t598 * t379, t601 * t369 - t751, t716, t768, t715, -pkin(5) * t337 + t601 * t300 - t598 * t316, -pkin(5) * t342 + t601 * t301 - t598 * t317, -pkin(5) * t327 + t601 * t293 - t598 * t299, -pkin(5) * t296 + t601 * t284 - t598 * t288, t698, t601 * t368 - t751, t601 * t357 - t598 * t377, t715, -t768, t716, -pkin(5) * t341 + t601 * t292 - t598 * t306, -pkin(5) * t328 + t601 * t287 - t598 * t294, -pkin(5) * t338 + t601 * t291 - t598 * t305, -pkin(5) * t289 + t601 * t282 - t598 * t283; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t560, 0, t561, 0, t528, t627, t496, pkin(5) * t496, t476, t446, t729, t475, -t721, t516, t601 * t443 + t598 * t452 - t753, t601 * t444 + t598 * t453 - t765, t598 * t440 + t486, pkin(5) * t401 - (-pkin(1) * t601 - pkin(6) * t598) * t440, t516, -t729, t721, t476, t446, t475, t598 * t370 + t601 * t550 + t486, t598 * t371 + t601 * t391 + t753, t598 * t372 + t601 * t385 + t765, pkin(5) * t351 + t598 * t335 + t601 * t346, t703, t598 * t356 + t601 * t379, t598 * t369 + t749, t718, t769, t717, pkin(5) * t339 + t598 * t300 + t601 * t316, pkin(5) * t344 + t598 * t301 + t601 * t317, pkin(5) * t329 + t598 * t293 + t601 * t299, pkin(5) * t297 + t598 * t284 + t601 * t288, t703, t598 * t368 + t749, t598 * t357 + t601 * t377, t717, -t769, t718, pkin(5) * t343 + t598 * t292 + t601 * t306, pkin(5) * t330 + t598 * t287 + t601 * t294, pkin(5) * t340 + t598 * t291 + t601 * t305, pkin(5) * t290 + t598 * t282 + t601 * t283; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t568, t569, 0, 0, t499, t483, t619, t498, t706, 0, pkin(1) * t557 + t664 - t723, -pkin(1) * t554 - t682 - t744, t441 + t657, pkin(1) * t532 + pkin(6) * t441, 0, -t619, -t706, t499, t483, t498, t600 * t427 + t597 * t428 + t657, t723 + t600 * t403 + (-pkin(1) - t689) * t557, t744 + t597 * t402 + (pkin(1) + t694) * t554, pkin(6) * t382 + (pkin(1) - t623) * t420, t626, t600 * t375 + t740, t597 * t438 + t750, t699, -t756, t702, pkin(6) * t360 + t597 * t321 + t600 * t323 - t747, pkin(6) * t365 + t597 * t322 + t600 * t325 - t767, -pkin(1) * t380 + pkin(6) * t354 + t600 * t307 + t597 * t336, -pkin(1) * t319 + pkin(6) * t315 + t600 * t298 + t597 * t302, t626, t597 * t434 + t750, t600 * t373 - t740, t702, t756, t699, pkin(6) * t364 + t600 * t309 + t597 * t313 - t747, -pkin(1) * t378 + pkin(6) * t355 + t600 * t295 + t597 * t320, pkin(6) * t361 + t600 * t308 + t597 * t312 + t767, -pkin(1) * t311 + pkin(6) * t304 + t600 * t285 + t597 * t286;];
tauB_reg = t1;
