% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:43:04
% EndTime: 2019-03-09 09:43:20
% DurationCPUTime: 11.58s
% Computational Cost: add. (6382->537), mult. (16628->741), div. (0->0), fcn. (12968->10), ass. (0->235)
t549 = sin(pkin(11));
t551 = cos(pkin(11));
t550 = sin(pkin(6));
t556 = sin(qJ(2));
t664 = t550 * t556;
t676 = cos(qJ(5));
t605 = t676 * t664;
t586 = qJD(2) * t605;
t573 = qJD(1) * t586;
t555 = sin(qJ(5));
t637 = qJD(1) * qJD(2);
t622 = t550 * t637;
t603 = t556 * t622;
t585 = t555 * t603;
t565 = t549 * t585 - t551 * t573;
t552 = cos(pkin(6));
t648 = qJD(1) * t552;
t535 = qJD(2) + t648;
t558 = cos(qJ(2));
t649 = qJD(1) * t550;
t628 = t558 * t649;
t480 = t535 * t551 - t549 * t628;
t654 = t549 * t535 + t551 * t628;
t688 = t480 * t676 - t555 * t654;
t396 = qJD(5) * t688 + t565;
t557 = cos(qJ(6));
t390 = t557 * t396;
t554 = sin(qJ(6));
t474 = t676 * t654;
t418 = t480 * t555 + t474;
t414 = qJD(6) + t418;
t692 = t414 ^ 2;
t693 = -t554 * t692 + t390;
t647 = qJD(1) * t556;
t629 = t550 * t647;
t519 = qJD(5) + t629;
t691 = t418 * t519;
t624 = qJD(5) * t676;
t690 = qJD(1) * t605 + t624;
t643 = qJD(5) * t555;
t689 = -t555 * t629 - t643;
t636 = pkin(1) * t648;
t487 = pkin(8) * t629 - t558 * t636;
t638 = qJD(3) + t487;
t656 = -t549 * t690 + t551 * t689;
t655 = t549 * t689 + t551 * t690;
t639 = pkin(3) * t629 + t638;
t547 = t556 ^ 2;
t685 = MDP(5) * (-t558 ^ 2 + t547);
t537 = pkin(8) * t664;
t631 = -pkin(1) * t558 - pkin(2);
t451 = pkin(3) * t664 + t537 + (-qJ(4) + t631) * t552;
t553 = -pkin(2) - qJ(4);
t621 = -qJ(3) * t556 - pkin(1);
t463 = (t553 * t558 + t621) * t550;
t407 = t551 * t451 - t463 * t549;
t663 = t550 * t558;
t491 = -t549 * t663 + t551 * t552;
t387 = pkin(4) * t664 - pkin(9) * t491 + t407;
t408 = t549 * t451 + t551 * t463;
t490 = t549 * t552 + t551 * t663;
t393 = -pkin(9) * t490 + t408;
t684 = t555 * t387 + t676 * t393;
t531 = pkin(2) * t629;
t595 = -qJ(3) * t558 + qJ(4) * t556;
t467 = t595 * t649 + t531;
t488 = pkin(8) * t628 + t556 * t636;
t472 = pkin(3) * t628 + t488;
t409 = -t467 * t549 + t551 * t472;
t665 = t549 * t556;
t568 = (pkin(4) * t558 - pkin(9) * t665) * t550;
t399 = qJD(1) * t568 + t409;
t410 = t551 * t467 + t549 * t472;
t662 = t551 * t556;
t611 = pkin(9) * t550 * t662;
t401 = qJD(1) * t611 + t410;
t575 = -t549 * t676 - t555 * t551;
t674 = -pkin(9) + t553;
t505 = t674 * t549;
t506 = t674 * t551;
t576 = -t555 * t505 + t506 * t676;
t683 = -qJD(4) * t575 - qJD(5) * t576 + t555 * t399 + t676 * t401;
t456 = t505 * t676 + t555 * t506;
t660 = t555 * t549;
t574 = t551 * t676 - t660;
t682 = -qJD(4) * t574 - qJD(5) * t456 - t399 * t676 + t555 * t401;
t675 = pkin(1) * t556;
t540 = t552 * t675;
t681 = pkin(8) * t663 + t540;
t630 = -pkin(4) * t551 - pkin(3);
t571 = t630 * t629;
t640 = -t571 + t638;
t521 = t558 * t622;
t428 = t535 * t553 + t639;
t449 = qJD(1) * t463;
t391 = t551 * t428 - t449 * t549;
t373 = pkin(4) * t629 - pkin(9) * t480 + t391;
t392 = t549 * t428 + t551 * t449;
t375 = -pkin(9) * t654 + t392;
t355 = t555 * t373 + t375 * t676;
t514 = pkin(2) * t603;
t644 = qJD(3) * t556;
t560 = (qJD(2) * t595 - qJD(4) * t558 - t644) * t550;
t422 = qJD(1) * t560 + t514;
t635 = pkin(1) * qJD(2) * t552;
t607 = qJD(1) * t635;
t481 = pkin(8) * t521 + t556 * t607;
t439 = pkin(3) * t521 - qJD(4) * t535 + t481;
t380 = -t422 * t549 + t551 * t439;
t563 = qJD(2) * t568;
t372 = qJD(1) * t563 + t380;
t381 = t551 * t422 + t549 * t439;
t597 = qJD(2) * t611;
t376 = qJD(1) * t597 + t381;
t561 = -qJD(5) * t355 + t676 * t372 - t555 * t376;
t347 = -pkin(5) * t521 - t561;
t680 = t414 * (pkin(5) * t688 + pkin(10) * t414) + t347;
t522 = t535 * qJ(3);
t442 = t522 + qJD(4) + t472;
t645 = qJD(2) * t558;
t650 = qJ(3) * qJD(2);
t679 = t556 * (qJD(4) - t442 + t650) - t553 * t645;
t646 = qJD(2) * t556;
t627 = t550 * t646;
t526 = pkin(2) * t627;
t433 = t526 + t560;
t677 = pkin(3) + pkin(8);
t453 = -qJD(4) * t552 + (t663 * t677 + t540) * qJD(2);
t397 = -t433 * t549 + t551 * t453;
t379 = t397 + t563;
t398 = t551 * t433 + t549 * t453;
t386 = t597 + t398;
t678 = -qJD(5) * t684 + t379 * t676 - t555 * t386;
t395 = -qJD(5) * t474 - t480 * t643 + t549 * t573 + t551 * t585;
t641 = qJD(6) * t557;
t632 = t557 * t395 + t519 * t641 + t554 * t521;
t642 = qJD(6) * t554;
t361 = -t642 * t688 + t632;
t673 = t361 * t554;
t672 = t391 * t558;
t671 = t392 * t558;
t668 = t688 * t554;
t404 = -t557 * t519 + t668;
t670 = t404 * t414;
t406 = t519 * t554 + t557 * t688;
t669 = t406 * t414;
t667 = t574 * t557;
t546 = t550 ^ 2;
t666 = t546 * qJD(1) ^ 2;
t661 = t554 * t396;
t541 = t549 * pkin(4) + qJ(3);
t657 = pkin(5) * t628 - t682;
t653 = pkin(8) * t603 - t558 * t607;
t530 = t558 * t635;
t543 = t552 * qJD(3);
t652 = t530 + t543;
t634 = t547 * t666;
t633 = t558 * t666;
t461 = -t535 * qJD(3) + t653;
t482 = -t552 * qJ(3) - t681;
t626 = t550 * t645;
t623 = t546 * t637;
t619 = t655 * t519;
t570 = t555 * t372 + t373 * t624 - t375 * t643 + t676 * t376;
t346 = pkin(10) * t521 + t570;
t416 = qJD(2) * t571 - t461;
t357 = pkin(5) * t396 - pkin(10) * t395 + t416;
t618 = -t346 * t554 + t557 * t357;
t616 = t395 * t554 - t557 * t521;
t615 = -t554 * t656 - t557 * t628;
t614 = t554 * t628 - t557 * t656;
t612 = t414 * t557;
t609 = -qJD(6) * t575 + t535;
t606 = t556 * t633;
t466 = pkin(3) * t663 - t482;
t602 = -0.2e1 * pkin(1) * t623;
t601 = t519 * t656 + t574 * t521;
t600 = t488 * t535 - t481;
t444 = -pkin(5) * t575 - pkin(10) * t574 + t541;
t599 = pkin(10) * t628 - qJD(6) * t444 + t683;
t598 = -pkin(5) * t655 + pkin(10) * t656 + qJD(6) * t456 - t640;
t594 = t346 * t557 + t357 * t554;
t353 = t519 * pkin(10) + t355;
t411 = pkin(4) * t654 + t442;
t365 = t418 * pkin(5) - pkin(10) * t688 + t411;
t349 = t353 * t557 + t365 * t554;
t593 = t353 * t554 - t365 * t557;
t359 = pkin(10) * t664 + t684;
t429 = pkin(4) * t490 + t466;
t437 = -t555 * t490 + t491 * t676;
t577 = -t490 * t676 - t555 * t491;
t368 = -pkin(5) * t577 - pkin(10) * t437 + t429;
t592 = t359 * t557 + t368 * t554;
t591 = -t359 * t554 + t368 * t557;
t590 = t380 * t551 + t381 * t549;
t589 = -t391 * t549 + t392 * t551;
t489 = t681 * qJD(2);
t588 = t481 * t552 + t489 * t535;
t584 = -pkin(8) * t627 + t530;
t583 = -t535 * t628 + t521;
t582 = -t437 * t554 + t557 * t664;
t413 = t437 * t557 + t554 * t664;
t354 = t373 * t676 - t555 * t375;
t581 = t387 * t676 - t555 * t393;
t569 = t555 * t379 + t676 * t386 + t387 * t624 - t393 * t643;
t483 = (-pkin(2) * t558 + t621) * t550;
t567 = -t574 * t641 + t615;
t566 = -t574 * t642 - t614;
t564 = (-qJ(3) * t645 - t644) * t550;
t438 = -pkin(3) * t603 - t461;
t352 = -t519 * pkin(5) - t354;
t562 = -pkin(10) * t396 + (t352 + t354) * t414;
t432 = (-pkin(8) + t630) * t627 + t652;
t486 = -qJ(3) * t628 + t531;
t485 = t552 * t631 + t537;
t476 = -t543 - t584;
t475 = qJD(1) * t483;
t473 = t526 + t564;
t468 = -t522 - t488;
t462 = -pkin(2) * t535 + t638;
t458 = qJD(1) * t564 + t514;
t452 = -t627 * t677 + t652;
t450 = t475 * t629;
t403 = qJD(5) * t437 - t551 * t586 + t627 * t660;
t402 = qJD(5) * t577 - t575 * t627;
t367 = qJD(6) * t413 + t402 * t554 - t557 * t626;
t366 = qJD(6) * t582 + t402 * t557 + t554 * t626;
t362 = qJD(6) * t406 + t616;
t360 = pkin(5) * t403 - pkin(10) * t402 + t432;
t358 = -pkin(5) * t664 - t581;
t351 = -pkin(5) * t626 - t678;
t350 = pkin(10) * t626 + t569;
t345 = -qJD(6) * t349 + t618;
t344 = -qJD(6) * t593 + t594;
t1 = [(MDP(6) * t626 - MDP(7) * t627) * (t535 + t648) + (-t461 * t558 + t481 * t556 + (t462 * t558 + t468 * t556) * qJD(2) + (-t476 * t558 + t489 * t556 + (t482 * t556 + t485 * t558) * qJD(2)) * qJD(1)) * t550 * MDP(11) + (t438 * t491 + t452 * t480 + ((-qJD(1) * t398 - t381) * t556 + (t442 * t665 - t671 + (-t408 * t558 + t466 * t665) * qJD(1)) * qJD(2)) * t550) * MDP(16) + (t452 * t654 + t438 * t490 + ((qJD(1) * t397 + t380) * t556 + (-t442 * t662 + t672 + (t407 * t558 - t466 * t662) * qJD(1)) * qJD(2)) * t550) * MDP(15) + (-t535 * t584 + t552 * t653 + t558 * t602) * MDP(10) + (-t398 * t654 - t381 * t490 - t397 * t480 - t380 * t491 + ((-t407 * t549 + t408 * t551) * qJD(1) + t589) * t627) * MDP(17) + ((-t475 * t646 + t458 * t558 + (t473 * t558 - t483 * t646) * qJD(1)) * t550 + t588) * MDP(12) + (-t461 * t552 - t476 * t535 + (-t475 * t645 - t458 * t556 + (-t473 * t556 - t483 * t645) * qJD(1)) * t550) * MDP(13) + (t556 * t602 - t588) * MDP(9) + (t458 * t483 + t461 * t482 + t462 * t489 + t468 * t476 + t473 * t475 + t481 * t485) * MDP(14) + (t380 * t407 + t381 * t408 + t391 * t397 + t392 * t398 + t438 * t466 + t442 * t452) * MDP(18) + (t361 * t413 + t366 * t406) * MDP(26) + 0.2e1 * (t556 * t558 * MDP(4) - t685) * t623 + (t361 * t582 - t362 * t413 - t366 * t404 - t367 * t406) * MDP(27) + (-t403 * t519 + (-t396 * t556 + (qJD(1) * t577 - t418) * t645) * t550) * MDP(22) + (-(qJD(6) * t591 + t350 * t557 + t360 * t554) * t414 - t592 * t396 + t344 * t577 - t349 * t403 + t351 * t406 + t358 * t361 + t347 * t413 + t352 * t366) * MDP(32) + (-t361 * t577 + t366 * t414 + t396 * t413 + t403 * t406) * MDP(28) + (-t396 * t577 + t403 * t414) * MDP(30) + (t354 * t626 + t429 * t396 + t411 * t403 - t416 * t577 + t432 * t418 + t519 * t678 + t581 * t521 + t561 * t664) * MDP(24) + (t362 * t577 - t367 * t414 + t396 * t582 - t403 * t404) * MDP(29) + ((-qJD(6) * t592 - t350 * t554 + t360 * t557) * t414 + t591 * t396 - t345 * t577 - t593 * t403 + t351 * t404 + t358 * t362 - t347 * t582 + t352 * t367) * MDP(31) + (t519 * t550 + t546 * t647) * MDP(23) * t645 + (t402 * t519 + (t395 * t556 + (qJD(1) * t437 + t688) * t645) * t550) * MDP(21) + (t395 * t437 + t402 * t688) * MDP(19) + (-t569 * t519 + t432 * t688 + t429 * t395 + t416 * t437 + t411 * t402 + (-t570 * t556 + (-qJD(1) * t684 - t355) * t645) * t550) * MDP(25) + (t395 * t577 - t396 * t437 - t402 * t418 - t403 * t688) * MDP(20); (-qJD(2) + t535) * MDP(7) * t629 + (t438 * t549 + (-t409 * t556 - t551 * t679 - t672) * t649 + t639 * t654) * MDP(15) + (t666 * t675 + t600) * MDP(9) + (pkin(1) * t633 - t487 * t535 + t653) * MDP(10) + (t410 * t654 + t409 * t480 + (qJD(4) * t480 - t392 * t629 - t380) * t551 + (qJD(4) * t654 + t391 * t629 - t381) * t549) * MDP(17) + (t638 * t535 + (t475 * t558 + t486 * t556) * t649 - t461) * MDP(13) + (-pkin(2) * t481 - qJ(3) * t461 - t462 * t488 - t468 * t638 - t475 * t486) * MDP(14) + (qJ(3) * t438 - t391 * t409 - t392 * t410 + t590 * t553 + t639 * t442 + (-t391 * t551 - t392 * t549) * qJD(4)) * MDP(18) + (-t486 * t628 + t450 - t600) * MDP(12) - MDP(4) * t606 + t583 * MDP(6) + t666 * t685 + (t438 * t551 + t639 * t480 + (t410 * t556 + t549 * t679 + t671) * t649) * MDP(16) + (t615 * t406 + t614 * t404 - (t673 + t362 * t557 + (-t404 * t554 + t406 * t557) * qJD(6)) * t574) * MDP(27) + (t362 * t575 - t404 * t655 + t414 * t567 - t574 * t661) * MDP(29) + (-t361 * t575 + t390 * t574 + t406 * t655 + t414 * t566) * MDP(28) + ((t444 * t557 - t456 * t554) * t396 - t345 * t575 - t576 * t362 + t347 * t554 * t574 + (t554 * t599 - t557 * t598) * t414 + t657 * t404 - t655 * t593 - t567 * t352) * MDP(31) + (t541 * t396 - t416 * t575 + t682 * t519 + t640 * t418 + t655 * t411 + (qJD(2) * t576 - t354) * t628) * MDP(24) + (-(t444 * t554 + t456 * t557) * t396 + t344 * t575 - t576 * t361 + t347 * t667 + (t554 * t598 + t557 * t599) * t414 + t657 * t406 - t655 * t349 + t566 * t352) * MDP(32) + (-t396 * t575 + t414 * t655) * MDP(30) + (-t619 + (qJD(2) * t575 + t418) * t628) * MDP(22) + (t361 * t667 + t406 * t566) * MDP(26) - t519 * MDP(23) * t628 + ((-t468 - t488 - t650) * t556 + (-pkin(2) * qJD(2) - t462 + t638) * t558) * MDP(11) * t649 + (-t628 * t688 + t601) * MDP(21) + (t395 * t574 + t656 * t688) * MDP(19) + (t541 * t395 + t416 * t574 + t683 * t519 + t640 * t688 + t656 * t411 + (-qJD(2) * t456 + t355) * t628) * MDP(25) + (t395 * t575 - t396 * t574 - t418 * t656 - t655 * t688) * MDP(20); t583 * MDP(11) + MDP(12) * t606 + (-t535 ^ 2 - t634) * MDP(13) + (t468 * t535 + t450 + t481) * MDP(14) + (t521 * t551 - t535 * t654 - t549 * t634) * MDP(15) + (-t480 * t535 - t521 * t549 - t551 * t634) * MDP(16) + (t549 * t480 - t551 * t654) * MDP(17) * t629 + (-t442 * t535 + t589 * t629 + t590) * MDP(18) + (-t418 * t535 + t601) * MDP(24) + (t521 * t575 - t535 * t688 - t619) * MDP(25) + (t575 * t661 - t574 * t362 - t656 * t404 + (-t554 * t655 - t557 * t609) * t414) * MDP(31) + (t575 * t390 - t574 * t361 - t656 * t406 + (t554 * t609 - t557 * t655) * t414) * MDP(32); (-t480 ^ 2 - t654 ^ 2) * MDP(17) + (t391 * t480 + t392 * t654 + t438) * MDP(18) + (t519 * t688 + t396) * MDP(24) + (t395 - t691) * MDP(25) + (-t404 * t688 + t693) * MDP(31) + (-t406 * t688 - t557 * t692 - t661) * MDP(32) + ((-qJD(2) * t551 + t480) * MDP(15) + (qJD(2) * t549 - t654) * MDP(16)) * t629; -t418 ^ 2 * MDP(20) + (t395 + t691) * MDP(21) - t565 * MDP(22) + MDP(23) * t521 + (t355 * t519 + t561) * MDP(24) + (t354 * t519 + t411 * t418 - t570) * MDP(25) + (t406 * t612 + t673) * MDP(26) + ((t361 - t670) * t557 + (-t362 - t669) * t554) * MDP(27) + (t414 * t612 + t661) * MDP(28) + t693 * MDP(29) + (-pkin(5) * t362 - t355 * t404 + t562 * t554 - t557 * t680) * MDP(31) + (-pkin(5) * t361 - t355 * t406 + t554 * t680 + t562 * t557) * MDP(32) + (t418 * MDP(19) + (-qJD(5) + t519) * MDP(22) - t411 * MDP(24) - t406 * MDP(28) + t404 * MDP(29) - t414 * MDP(30) + t593 * MDP(31) + t349 * MDP(32) + MDP(20) * t688) * t688; t406 * t404 * MDP(26) + (-t404 ^ 2 + t406 ^ 2) * MDP(27) + (t632 + t670) * MDP(28) + (-t616 + t669) * MDP(29) + t396 * MDP(30) + (t349 * t414 - t352 * t406 + t618) * MDP(31) + (t352 * t404 - t414 * t593 - t594) * MDP(32) + (-MDP(28) * t668 - MDP(29) * t406 - MDP(31) * t349 + MDP(32) * t593) * qJD(6);];
tauc  = t1;
