% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:24:00
% EndTime: 2019-03-09 08:24:10
% DurationCPUTime: 9.21s
% Computational Cost: add. (3253->651), mult. (7253->777), div. (0->0), fcn. (4859->8), ass. (0->281)
t742 = pkin(4) + qJ(3);
t561 = cos(pkin(9));
t566 = cos(qJ(2));
t662 = qJD(1) * qJD(2);
t638 = t566 * t662;
t563 = sin(qJ(2));
t660 = qJDD(1) * t563;
t586 = t638 + t660;
t560 = sin(pkin(9));
t659 = qJDD(2) * t560;
t447 = t561 * t586 + t659;
t534 = pkin(7) * t660;
t459 = -qJDD(2) * pkin(2) + pkin(7) * t638 + qJDD(3) + t534;
t678 = qJD(1) * t563;
t645 = t561 * t678;
t676 = qJD(2) * t560;
t483 = t645 + t676;
t670 = qJD(4) * t483;
t540 = t561 * qJDD(2);
t446 = t560 * t586 - t540;
t718 = pkin(3) * t446;
t574 = t459 - t670 + t718;
t711 = -pkin(5) - qJ(4);
t647 = t560 * t678;
t665 = t561 * qJD(2);
t481 = t647 - t665;
t733 = t446 * qJ(5) + t481 * qJD(5);
t374 = t447 * t711 + t574 + t733;
t564 = sin(qJ(1));
t567 = cos(qJ(1));
t621 = g(1) * t567 + g(2) * t564;
t724 = t563 * t621;
t741 = -t374 + t724;
t694 = t563 * t567;
t695 = t563 * t564;
t740 = g(1) * t694 + g(2) * t695;
t739 = MDP(15) - MDP(20);
t738 = MDP(17) + MDP(19);
t737 = -2 * pkin(1);
t712 = pkin(3) + qJ(5);
t736 = t481 * t712;
t735 = t561 * t712;
t661 = qJD(1) * qJD(4);
t663 = qJ(4) * qJDD(1);
t734 = t566 * (t661 + t663);
t557 = g(3) * t566;
t634 = t459 + t557;
t669 = qJD(4) * t560;
t677 = qJD(1) * t566;
t537 = pkin(7) * t677;
t646 = t560 * t677;
t684 = pkin(3) * t646 + t537;
t732 = -qJD(5) * t561 - t669 - t684;
t549 = t566 * qJ(5);
t698 = t561 * t563;
t731 = pkin(4) * t698 + t549;
t544 = t566 * qJDD(1);
t545 = t566 * qJD(5);
t730 = qJ(5) * t544 + qJD(1) * t545;
t639 = t563 * t662;
t729 = -t639 + t544;
t675 = qJD(2) * t563;
t728 = -qJ(4) * t675 + qJD(4) * t566;
t725 = -MDP(16) + MDP(21);
t674 = qJD(2) * t566;
t643 = t560 * t674;
t683 = pkin(3) * t643 + pkin(7) * t674;
t723 = -t563 * (qJD(4) * t561 - qJD(5) * t560) + t683;
t477 = t481 ^ 2;
t721 = t483 ^ 2;
t720 = pkin(4) + pkin(8);
t569 = qJD(1) ^ 2;
t719 = pkin(1) * t569;
t717 = pkin(4) * t447;
t716 = pkin(7) * t483;
t715 = g(1) * t564;
t553 = t566 * pkin(2);
t710 = qJ(4) * t447;
t709 = qJ(5) * t560;
t565 = cos(qJ(6));
t562 = sin(qJ(6));
t704 = t481 * t562;
t422 = -t565 * t483 + t704;
t521 = -qJD(6) + t677;
t708 = t422 * t521;
t424 = t481 * t565 + t483 * t562;
t707 = t424 * t521;
t615 = pkin(2) * t563 - qJ(3) * t566;
t463 = qJD(2) * t615 - qJD(3) * t563;
t706 = t463 * t561;
t705 = t481 * t483;
t490 = t615 * qJD(1);
t703 = t490 * t561;
t559 = t566 ^ 2;
t702 = t559 * t569;
t701 = t560 * qJ(4);
t700 = t560 * t563;
t699 = t560 * t566;
t543 = t561 * qJ(3);
t697 = t561 * t565;
t696 = t561 * t566;
t547 = t563 * qJ(3);
t693 = t564 * t566;
t692 = t566 * t567;
t691 = t567 * t560;
t650 = -pkin(1) - t553;
t599 = t650 - t547;
t419 = qJD(1) * t463 + qJDD(1) * t599;
t454 = pkin(7) * t729 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t389 = t560 * t419 + t561 * t454;
t486 = t560 * t565 + t561 * t562;
t465 = t486 * qJD(6);
t587 = t566 * t486;
t690 = -qJD(1) * t587 + t465;
t487 = t560 * t562 - t697;
t466 = t487 * qJD(6);
t644 = t561 * t677;
t689 = -t562 * t646 + t565 * t644 + t466;
t472 = t560 * t490;
t531 = qJ(4) * t678;
t688 = t472 + t531;
t474 = t599 * qJD(1);
t503 = qJD(2) * qJ(3) + t537;
t427 = t560 * t474 + t561 * t503;
t593 = t561 * t711 + t709;
t577 = t566 * t593;
t687 = -qJD(1) * t577 + t732;
t611 = qJ(4) * t561 - t709;
t588 = t566 * t611;
t686 = -qJD(1) * t588 - t732;
t685 = t740 * t560;
t682 = pkin(3) * t700 + t563 * pkin(7);
t680 = t553 + t547;
t494 = -pkin(1) - t680;
t449 = pkin(7) * t696 + t560 * t494;
t681 = g(1) * t695 - g(2) * t694;
t497 = t742 * t560;
t498 = t561 * pkin(4) + t543;
t558 = t563 ^ 2;
t679 = t558 - t559;
t673 = qJD(3) * t481;
t672 = qJD(3) * t483;
t671 = qJD(3) * t561;
t667 = qJD(6) * t562;
t666 = qJD(6) * t565;
t426 = t474 * t561 - t560 * t503;
t412 = pkin(3) * t677 + qJD(4) - t426;
t584 = qJ(5) * t677 + t412;
t392 = pkin(4) * t483 + t584;
t664 = -qJD(4) - t392;
t656 = pkin(7) * t675;
t655 = t560 * t720;
t654 = -qJ(4) * t639 - t389;
t653 = t565 * t446 + t562 * t447 + t483 * t666;
t652 = qJD(3) * t646 + t561 * t740;
t651 = g(1) * t692 + g(2) * t693 + g(3) * t563;
t649 = -pkin(7) * t560 - pkin(3);
t648 = qJ(3) * t675;
t642 = t566 * t665;
t641 = t711 * t566;
t640 = qJ(3) * t544;
t637 = pkin(2) + t701;
t383 = t574 - t710;
t377 = -t383 - t733;
t636 = t377 - t557;
t635 = t383 + t557;
t388 = t419 * t561 - t560 * t454;
t596 = pkin(3) * t544 + qJDD(4) - t388;
t571 = -t639 * t712 + t596 + t730;
t375 = t447 * t720 + t571;
t626 = qJDD(5) - t654;
t376 = pkin(5) * t639 - t720 * t446 + (qJDD(1) * t711 - t661) * t566 + t626;
t633 = -t562 * t375 + t565 * t376;
t632 = t562 * t446 - t565 * t447;
t630 = pkin(4) * t644 - t703;
t513 = pkin(7) * t699;
t448 = t494 * t561 - t513;
t628 = pkin(3) * t696 + qJ(4) * t699 + t680;
t627 = t567 * pkin(1) + pkin(2) * t692 + t564 * pkin(7) + qJ(3) * t694;
t625 = t649 * t563;
t624 = -pkin(2) - t735;
t468 = t560 * t693 + t561 * t567;
t470 = -t564 * t561 + t566 * t691;
t623 = g(1) * t468 - g(2) * t470;
t469 = t561 * t693 - t691;
t471 = t560 * t564 + t561 * t692;
t622 = g(1) * t469 - g(2) * t471;
t620 = -g(2) * t567 + t715;
t619 = qJD(2) * pkin(2) - pkin(7) * t678 - qJD(3);
t554 = t567 * pkin(7);
t618 = -t469 * pkin(3) - qJ(4) * t468 + t554;
t435 = qJ(4) * t566 - t449;
t617 = pkin(4) * t642 + t545 - t706;
t552 = t566 * pkin(3);
t439 = -t448 + t552;
t452 = t560 * t463;
t616 = t452 - t728;
t614 = pkin(7) * t481 - t560 * t619;
t612 = t560 * t660 - t540;
t610 = t565 * t375 + t562 * t376;
t385 = t483 * t720 + t584;
t386 = qJD(1) * t641 - t481 * t720 + qJD(5) + t427;
t372 = t385 * t565 + t386 * t562;
t609 = t385 * t562 - t386 * t565;
t402 = t513 + t552 + (pkin(8) * t563 - t494) * t561 + t731;
t403 = -t563 * t655 + t449 + t641;
t608 = t402 * t565 + t403 * t562;
t607 = -t402 * t562 + t403 * t565;
t606 = t468 * t565 + t469 * t562;
t605 = t468 * t562 - t469 * t565;
t604 = -qJ(3) * t446 - t673;
t603 = qJ(3) * t447 + t672;
t568 = qJD(2) ^ 2;
t601 = qJDD(2) * t566 - t563 * t568;
t441 = -pkin(7) * t645 + t472;
t433 = -t561 * t656 + t452;
t600 = (-qJ(5) + t649) * t563;
t598 = pkin(3) * t561 + t637;
t417 = qJ(4) * t677 - t427;
t485 = qJDD(6) - t729;
t595 = -t485 * t562 + t521 * t666;
t594 = t485 * t565 + t521 * t667;
t476 = pkin(8) * t561 + t498;
t576 = pkin(8) * t696 + t600;
t592 = -qJD(1) * t576 + qJD(3) * t560 + qJD(6) * t476 - t630;
t475 = pkin(8) * t560 + t497;
t575 = (-pkin(7) * t561 + pkin(5)) * t563 - t566 * t655;
t591 = -qJD(1) * t575 - qJD(6) * t475 + t671 - t688;
t590 = -pkin(4) * t699 - pkin(7) * t698;
t589 = -pkin(4) * t446 + t626;
t585 = t561 * t660 + t659;
t380 = -t481 * t667 + t653;
t583 = t561 * t640 - t685;
t581 = qJ(4) * t483 + t619;
t580 = t471 * pkin(3) + qJ(4) * t470 + t627;
t579 = t599 * t715;
t578 = t560 * t640 + t652;
t381 = qJD(6) * t424 + t632;
t573 = t634 - t724;
t572 = -g(1) * t470 - g(2) * t468 - g(3) * t700 + t596;
t570 = t573 - t710 + t718;
t519 = qJ(3) * t692;
t514 = qJ(3) * t693;
t467 = t637 + t735;
t458 = t486 * t563;
t457 = t562 * t700 - t563 * t697;
t455 = -qJ(4) * t698 + t682;
t453 = t560 * t711 + t624;
t445 = -qJ(4) * t644 + t684;
t440 = pkin(7) * t647 + t703;
t434 = t563 * t611 - t682;
t432 = t560 * t656 + t706;
t431 = qJD(1) * t625 - t703;
t430 = -t441 - t531;
t429 = (-qJ(4) * t674 - qJD(4) * t563) * t561 + t683;
t425 = -pkin(4) * t700 - t435;
t421 = t563 * t593 + t682;
t420 = qJD(2) * t625 - t706;
t418 = t439 + t731;
t416 = t470 * t565 + t471 * t562;
t415 = -t470 * t562 + t471 * t565;
t414 = qJD(1) * t590 + t688;
t411 = -t433 + t728;
t407 = pkin(3) * t481 - t581;
t405 = qJD(2) * t587 - t466 * t563;
t404 = t465 * t563 + t562 * t643 - t565 * t642;
t401 = qJD(1) * t600 + t630;
t400 = qJD(2) * t588 - t723;
t399 = qJD(2) * t590 + t616;
t396 = qJD(2) * t600 + t617;
t395 = -pkin(4) * t481 + qJD(5) - t417;
t394 = t581 - t736;
t393 = qJD(2) * t577 + t723;
t391 = qJD(2) * t575 + t616;
t390 = qJD(2) * t576 + t617;
t387 = t483 * t711 - t619 + t736;
t384 = -pkin(3) * t639 + t596;
t382 = t654 + t734;
t379 = t589 - t734;
t378 = t571 + t717;
t1 = [(qJDD(2) * t563 + t566 * t568) * MDP(6) + ((-pkin(7) * qJDD(2) + t662 * t737) * t563 + (0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t568 + t620) * t566) * MDP(9) + (-pkin(7) * t601 + t586 * t737 - t681) * MDP(10) + (-t380 * t457 - t381 * t458 - t404 * t424 - t405 * t422) * MDP(24) + (t380 * t458 + t405 * t424) * MDP(23) + (qJDD(1) * t558 + 0.2e1 * t563 * t638) * MDP(4) + (-t429 * t483 - t447 * t455 + (-t383 * t561 + (-qJD(1) * t435 - t417) * qJD(2)) * t563 + (qJD(1) * t411 + qJDD(1) * t435 - t407 * t665 + t382) * t566 + t623) * MDP(17) + (t400 * t483 + t434 * t447 + (t377 * t561 + (qJD(1) * t425 + t395) * qJD(2)) * t563 + (-qJD(1) * t399 - qJDD(1) * t425 + t394 * t665 - t379) * t566 + t623) * MDP(19) + 0.2e1 * (t544 * t563 - t662 * t679) * MDP(5) + (-t432 * t483 - t433 * t481 - t446 * t449 - t447 * t448 + (-t388 * t561 - t389 * t560) * t563 + (-t426 * t561 - t427 * t560) * t674 + t681) * MDP(13) + (t411 * t481 + t420 * t483 + t435 * t446 + t439 * t447 + (t382 * t560 + t384 * t561) * t563 + (t412 * t561 + t417 * t560) * t674 + t681) * MDP(15) + (-t396 * t483 + t399 * t481 - t418 * t447 + t425 * t446 + (-t378 * t561 + t379 * t560) * t563 + (-t392 * t561 + t395 * t560) * t674 - t681) * MDP(20) + (-t485 * t566 - t521 * t675) * MDP(27) + (t381 * t566 + t404 * t521 - t422 * t675 - t457 * t485) * MDP(26) + (-t380 * t566 - t405 * t521 + t424 * t675 + t458 * t485) * MDP(25) + (-t429 * t481 - t446 * t455 + (-t383 * t560 + (qJD(1) * t439 + t412) * qJD(2)) * t563 + (-qJD(1) * t420 - qJDD(1) * t439 - t407 * t676 - t384) * t566 - t622) * MDP(16) + (-t400 * t481 - t434 * t446 + (-t377 * t560 + (-qJD(1) * t418 - t392) * qJD(2)) * t563 + (qJD(1) * t396 + qJDD(1) * t418 - t394 * t676 + t378) * t566 + t622) * MDP(21) + ((t390 * t565 + t391 * t562) * t521 - t608 * t485 + t610 * t566 - t372 * t675 + t393 * t424 + t421 * t380 + t374 * t458 + t387 * t405 - g(1) * t605 - g(2) * t415 + (t521 * t607 - t566 * t609) * qJD(6)) * MDP(29) + ((pkin(7) * t446 + t459 * t560 + (qJD(1) * t448 + t426) * qJD(2)) * t563 + (-qJD(1) * t432 + qJD(2) * t614 - qJDD(1) * t448 - t388) * t566 + t622) * MDP(11) + (-g(1) * t618 - g(2) * t580 + t382 * t435 + t383 * t455 + t384 * t439 + t407 * t429 + t417 * t411 + t412 * t420 - t579) * MDP(18) + t620 * MDP(2) + t621 * MDP(3) + t601 * MDP(7) + (-(-t390 * t562 + t391 * t565) * t521 + t607 * t485 - t633 * t566 - t609 * t675 + t393 * t422 + t421 * t381 + t374 * t457 + t387 * t404 + g(1) * t606 - g(2) * t416 + (t372 * t566 + t521 * t608) * qJD(6)) * MDP(28) + ((pkin(7) * t447 + t459 * t561 + (-qJD(1) * t449 - t427) * qJD(2)) * t563 + (qJD(1) * t433 + qJDD(1) * t449 + t389 + (-t561 * t619 + t716) * qJD(2)) * t566 - t623) * MDP(12) + (t389 * t449 + t427 * t433 + t388 * t448 + t426 * t432 - g(1) * t554 - g(2) * t627 - t579 + (t459 * t563 - t619 * t674) * pkin(7)) * MDP(14) + qJDD(1) * MDP(1) + (t378 * t418 + t392 * t396 + t377 * t434 + t394 * t400 + t379 * t425 + t395 * t399 - g(1) * (-qJ(5) * t469 + t618) - g(2) * (pkin(4) * t694 + qJ(5) * t471 + t580) - (-t563 * t742 + t650) * t715) * MDP(22); (-MDP(4) * t563 * t566 + MDP(5) * t679) * t569 + (t378 * t497 + t377 * t467 + t379 * t498 - t392 * t401 - t395 * t414 - g(1) * (pkin(4) * t692 + t519) - g(2) * (pkin(4) * t693 + t514) - g(3) * (t549 * t561 + t628) + t686 * t394 + (t392 * t560 + t395 * t561) * qJD(3) + (-g(3) * pkin(4) + t621 * (-t624 + t701)) * t563) * MDP(22) + ((-pkin(7) * qJDD(1) + t719) * t566 + t651) * MDP(10) + (-t534 - t557 + (t621 + t719) * t563) * MDP(9) + (-t498 * t544 + t447 * t467 + t636 * t560 + t686 * t483 + ((qJD(2) * t498 - t395) * t563 + (t414 + (-qJD(3) - t394) * t561) * t566) * qJD(1) + t685) * MDP(19) + (t497 * t544 - t446 * t467 + t636 * t561 - t686 * t481 + ((t394 * t560 - t401) * t566 + (-qJD(2) * t497 + t392) * t563) * qJD(1) + t652) * MDP(21) + (t422 * t678 + t485 * t486 + t521 * t689) * MDP(26) + (t380 * t487 + t424 * t690) * MDP(23) + (t380 * t486 - t381 * t487 - t422 * t690 - t424 * t689) * MDP(24) + (-t424 * t678 + t485 * t487 - t521 * t690) * MDP(25) + (t440 * t483 + t441 * t481 + (t426 * t677 + t389 + t604) * t561 + (t427 * t677 - t388 + t603) * t560 - t651) * MDP(13) + (t401 * t483 - t414 * t481 + t446 * t498 - t447 * t497 + (t392 * t677 - t379 + t673) * t561 + (-t395 * t677 - t378 - t672) * t560 + t651) * MDP(20) + (-pkin(2) * t446 - t634 * t561 + ((-qJ(3) * t676 - t426) * t563 + (t440 - t614) * t566) * qJD(1) + t578) * MDP(11) + (-t430 * t481 - t431 * t483 + (-t412 * t677 - t382 + t604) * t561 + (-t417 * t677 + t384 + t603) * t560 - t651) * MDP(15) + t521 * MDP(27) * t678 + (-t382 * t543 - t407 * t445 - t412 * t431 - g(1) * t519 - g(2) * t514 - g(3) * t628 + (-t430 - t671) * t417 + (qJ(3) * t384 + qJD(3) * t412 - qJD(4) * t407) * t560 + (-t383 + t724) * t598) * MDP(18) + (t446 * t598 + t635 * t561 + (t445 + t669) * t481 + (-t412 * t563 + t431 * t566 + (t407 * t566 + t648) * t560) * qJD(1) - t578) * MDP(16) + (t445 * t483 + t447 * t598 + (-t635 + t670) * t560 + (t417 * t563 - t430 * t566 + (t648 + (-qJD(3) + t407) * t566) * t561) * qJD(1) - t583) * MDP(17) + (-pkin(2) * t447 + t634 * t560 + ((-qJ(3) * t665 + t427) * t563 + (-t716 - t441 + (qJD(3) + t619) * t561) * t566) * qJD(1) + t583) * MDP(12) + (-t459 * pkin(2) - t427 * t441 - t426 * t440 + t619 * t537 - g(1) * (-pkin(2) * t694 + t519) - g(2) * (-pkin(2) * t695 + t514) - g(3) * t680 + (-t426 * t560 + t427 * t561) * qJD(3) + (-t388 * t560 + t389 * t561) * qJ(3)) * MDP(14) + qJDD(2) * MDP(8) + (-(t475 * t565 + t476 * t562) * t485 + t453 * t380 + (t562 * t591 + t565 * t592) * t521 + t687 * t424 + t690 * t387 + t372 * t678 + (t557 - t741) * t487) * MDP(29) + ((-t475 * t562 + t476 * t565) * t485 + t453 * t381 - g(3) * t587 + (t562 * t592 - t565 * t591) * t521 + t687 * t422 + t689 * t387 + t609 * t678 + t741 * t486) * MDP(28) + MDP(7) * t544 + MDP(6) * t660; (t426 * t483 + t427 * t481 + t573) * MDP(14) + (-t417 * t481 + (-qJD(4) - t412) * t483 + t570) * MDP(18) + (t395 * t481 + t483 * t664 + t570 + t733) * MDP(22) + (t381 - t707) * MDP(28) + (t380 + t708) * MDP(29) + (MDP(11) + t725) * ((-t483 + t676) * t677 + t612) + (MDP(12) - t738) * ((t481 + t665) * t677 + t585) + (-MDP(13) - t739) * (t721 + t477); (t407 * t483 + t572) * MDP(18) + (-t394 * t483 + t572 + t717 + t730) * MDP(22) + (t422 * t483 + t595) * MDP(28) + (t424 * t483 - t594) * MDP(29) + t739 * ((-t481 + t665) * t677 + t585) + ((-t417 * MDP(18) + t395 * MDP(22) + (-t565 * MDP(28) + t562 * MDP(29)) * t521) * t566 + (-pkin(3) * MDP(18) - MDP(22) * t712 - t725) * t675) * qJD(1) + t725 * (t544 + t705) + t738 * (-t721 - t702); (-t729 + t705) * MDP(19) + t612 * MDP(20) + (-t477 - t702) * MDP(21) + (-g(1) * t471 - g(2) * t469 - g(3) * t698 + t394 * t481 + t589) * MDP(22) + (-t422 * t481 + t594) * MDP(28) + (-t424 * t481 + t595) * MDP(29) + (-MDP(22) * t663 + ((t483 + t676) * MDP(20) + t664 * MDP(22) + (-t562 * MDP(28) - t565 * MDP(29)) * t521) * qJD(1)) * t566; t424 * t422 * MDP(23) + (-t422 ^ 2 + t424 ^ 2) * MDP(24) + (t653 - t708) * MDP(25) + (-t632 - t707) * MDP(26) + t485 * MDP(27) + (-g(1) * t415 + g(2) * t605 + g(3) * t457 - t372 * t521 - t387 * t424 + t633) * MDP(28) + (g(1) * t416 + g(2) * t606 + g(3) * t458 + t387 * t422 + t521 * t609 - t610) * MDP(29) + (-MDP(25) * t704 - MDP(26) * t424 - MDP(28) * t372 + MDP(29) * t609) * qJD(6);];
tau  = t1;
